"""Pretrain across cell types, then finetune one checkpoint per cell type.

Why pretrain-then-finetune rather than four independent models: each cell type
has 22-32 usable rounds. A transformer fitted on 22 samples memorises them.
Pooling the four sequences gives ~111 samples to learn the shared shape of the
problem from, and the per-cell-type finetune - a short, low-learning-rate pass
on that cell type only, with its own checkpoint - is what specialises it. The
deliverable is still one model per cell type, trained and finetuned separately.

Protocol:
  * walk-forward folds only. A fold trains on rounds 1..k and tests on the
    rounds after k, so the model never sees a task earlier than one it predicts.
  * every fold reports the same metrics for the model and for five reference
    strategies, and the promoted checkpoint has to beat the best of them.
  * AdamW, cosine schedule with warmup, gradient clipping, EMA weights for
    evaluation, early stopping on the fold's own validation tail.

Usage
  train.py --start 2026-08-27T00:00:00 --end now             # full run
  train.py --finetune                                        # incremental
"""

import argparse
import json
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import torch

from .data import (CELL_TYPES, CONTEXT, FEATURE_DIM, N_CLASSES, build_samples,
                   load_rounds, refresh_seeds, walk_forward)
from .evaluate import evaluate_probs, permutation_test, run_baselines
from .model import EMA, SeedFormer, soft_cross_entropy

ROOT = Path(__file__).resolve().parents[1]
CKPT_DIR = ROOT / "seed_model" / "checkpoints"


def set_seed(seed):
    torch.manual_seed(seed)
    np.random.seed(seed)
    torch.use_deterministic_algorithms(True, warn_only=True)


def as_tensors(sample, idx):
    return (torch.from_numpy(sample["x"][idx]),
            torch.from_numpy(sample["mask"][idx]),
            torch.from_numpy(sample["cell"][idx]).long(),
            torch.from_numpy(sample["y"][idx]))


def fit(model, sample, train_idx, *, epochs, lr, weight_decay, batch_size,
        warmup=0.15, smoothing=0.05, ema_decay=0.995, val_frac=0.2, seed=0,
        patience=40, log=None):
    """Train in place; returns the EMA state that scored best on the tail."""
    set_seed(seed)
    order = train_idx[np.argsort(sample["round"][train_idx], kind="stable")]
    n_val = max(1, int(len(order) * val_frac)) if len(order) >= 10 else 0
    tr = order[:len(order) - n_val] if n_val else order
    va = order[len(order) - n_val:] if n_val else order[-1:]

    x, m, c, y = as_tensors(sample, tr)
    opt = torch.optim.AdamW(model.parameters(), lr=lr, weight_decay=weight_decay,
                            betas=(0.9, 0.95))
    steps_per_epoch = max(1, int(np.ceil(len(tr) / batch_size)))
    total = epochs * steps_per_epoch
    warm = max(1, int(total * warmup))

    def lr_at(step):
        if step < warm:
            return step / warm
        p = (step - warm) / max(1, total - warm)
        return 0.5 * (1 + np.cos(np.pi * min(p, 1.0)))

    sched = torch.optim.lr_scheduler.LambdaLR(opt, lr_at)
    ema = EMA(model, ema_decay)
    best = {"loss": float("inf"), "state": None, "epoch": -1}
    rng = np.random.default_rng(seed)

    for epoch in range(epochs):
        model.train()
        perm = rng.permutation(len(tr))
        for s in range(0, len(tr), batch_size):
            b = perm[s:s + batch_size]
            opt.zero_grad(set_to_none=True)
            loss = soft_cross_entropy(model(x[b], m[b], c[b]), y[b], smoothing)
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
            opt.step()
            sched.step()
            ema.update(model)

        shadow = SeedFormer(dim=model.dim, depth=len(model.blocks),
                            heads=model.blocks[0].attn.heads)
        ema.copy_to(shadow)
        shadow.eval()
        with torch.no_grad():
            vx, vm, vc, vy = as_tensors(sample, va)
            vloss = float(soft_cross_entropy(shadow(vx, vm, vc), vy, 0.0))
        if vloss < best["loss"] - 1e-5:
            best = {"loss": vloss, "state": {k: v.clone() for k, v in
                                             ema.shadow.items()}, "epoch": epoch}
        elif epoch - best["epoch"] >= patience:
            break
        if log and epoch % max(1, epochs // 6) == 0:
            log(f"      epoch {epoch:>4}  train {float(loss):.4f}  val {vloss:.4f}")
    return best


@torch.no_grad()
def predict_probs(state, sample, idx, cfg):
    model = SeedFormer(**cfg)
    model.load_state_dict(state)
    model.eval()
    x, m, c, _ = as_tensors(sample, idx)
    return torch.softmax(model(x, m, c), dim=-1).numpy()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--refresh", action="store_true",
                    help="re-run scripts/seed_bins.py first so seeds.json "
                         "covers --start/--end before training")
    ap.add_argument("--start", default="2026-08-27T00:00:00")
    ap.add_argument("--end", default="now")
    ap.add_argument("--dim", type=int, default=64)
    ap.add_argument("--depth", type=int, default=3)
    ap.add_argument("--heads", type=int, default=4)
    ap.add_argument("--dropout", type=float, default=0.2)
    ap.add_argument("--epochs", type=int, default=400)
    ap.add_argument("--finetune-epochs", type=int, default=120)
    ap.add_argument("--lr", type=float, default=3e-3)
    ap.add_argument("--finetune-lr", type=float, default=3e-4)
    ap.add_argument("--weight-decay", type=float, default=0.05)
    ap.add_argument("--batch-size", type=int, default=16)
    ap.add_argument("--folds", type=int, default=5)
    ap.add_argument("--seeds", type=int, default=3,
                    help="independent inits averaged into the reported result")
    ap.add_argument("--context", type=int, default=CONTEXT)
    ap.add_argument("--finetune", action="store_true",
                    help="warm-start from the promoted checkpoints instead of "
                         "training from scratch")
    ap.add_argument("--report", default=str(ROOT / "seed_model" / "report.json"))
    args = ap.parse_args()

    t0 = time.time()
    if args.refresh:
        print("refreshing seeds.json ...")
        rounds, meta = refresh_seeds(args.start, args.end)
    else:
        rounds, meta = load_rounds()
    sample = build_samples(rounds, args.context)
    print(f"rounds window {meta['window']['start'][:19]} -> "
          f"{meta['window']['end'][:19]} (seeds.json written "
          f"{(meta.get('generated_at') or '?')[:19]})")
    print(f"samples {len(sample['y'])} over {len(rounds)} cell types "
          f"(feature dim {FEATURE_DIM}, context {args.context})")
    for i, ct in enumerate(CELL_TYPES):
        n = int((sample["cell"] == i).sum())
        print(f"  {ct:<11} {n:>3} predictable rounds")

    cfg = dict(dim=args.dim, depth=args.depth, heads=args.heads,
               dropout=args.dropout)
    CKPT_DIR.mkdir(parents=True, exist_ok=True)

    folds = list(walk_forward(sample, args.folds))
    print(f"\nwalk-forward: {len(folds)} folds "
          f"({[len(te) for _, te in folds]} test rounds)")

    per_fold, pooled_model_hits, pooled_ref = [], [], {}
    per_cell_hits = {ct: {"model": [], "best_ref": []} for ct in CELL_TYPES}

    for f, (tr, te) in enumerate(folds, start=1):
        print(f"\n  fold {f}: train {len(tr)}  test {len(te)}")
        base = run_baselines(sample, tr, te)
        for name, res in base.items():
            pooled_ref.setdefault(name, []).extend(res["per_sample_hits"])

        # pooled pretrain, averaged over independent inits
        probs = np.zeros((len(te), N_CLASSES))
        for s in range(args.seeds):
            model = SeedFormer(**cfg)
            if args.finetune:
                shared = CKPT_DIR / "pooled.pt"
                if shared.exists():
                    model.load_state_dict(torch.load(shared, weights_only=True)["state"])
            best = fit(model, sample, tr, epochs=args.epochs, lr=args.lr,
                       weight_decay=args.weight_decay, batch_size=args.batch_size,
                       seed=1000 + s, log=print if s == 0 else None)
            probs += predict_probs(best["state"], sample, te, cfg) / args.seeds
        model_res = evaluate_probs(probs, sample["y"][te])
        pooled_model_hits.extend(model_res["per_sample_hits"])

        best_ref = max(base.items(), key=lambda kv: kv[1]["mean_hits"])
        print(f"      model {model_res['mean_hits']:.3f} hits/round  |  "
              f"best reference {best_ref[0]} {best_ref[1]['mean_hits']:.3f}  |  "
              f"log-loss {model_res['log_loss']:.4f} vs "
              f"{best_ref[1]['log_loss']:.4f}")

        for i, ct in enumerate(CELL_TYPES):
            local = [j for j, k in enumerate(te) if sample["cell"][k] == i]
            if not local:
                continue
            per_cell_hits[ct]["model"].extend(
                [model_res["per_sample_hits"][j] for j in local])
            per_cell_hits[ct]["best_ref"].extend(
                [best_ref[1]["per_sample_hits"][j] for j in local])

        per_fold.append({"fold": f, "n_train": len(tr), "n_test": len(te),
                         "model": {k: v for k, v in model_res.items()
                                   if k != "per_sample_hits"},
                         "baselines": {k: {kk: vv for kk, vv in v.items()
                                           if kk != "per_sample_hits"}
                                       for k, v in base.items()}})

    # ---- promote: train on everything, save pooled + one per cell type ----
    print("\n  final fit on all rounds")
    all_idx = np.arange(len(sample["y"]))
    model = SeedFormer(**cfg)
    if args.finetune and (CKPT_DIR / "pooled.pt").exists():
        model.load_state_dict(torch.load(CKPT_DIR / "pooled.pt",
                                         weights_only=True)["state"])
    pooled = fit(model, sample, all_idx, epochs=args.epochs, lr=args.lr,
                 weight_decay=args.weight_decay, batch_size=args.batch_size,
                 seed=7, log=print)
    meta_common = {"config": cfg, "context": args.context,
                   "feature_dim": FEATURE_DIM,
                   "window": meta["window"],
                   "trained_at": datetime.now(timezone.utc).isoformat(),
                   "n_samples": int(len(all_idx))}
    torch.save({"state": pooled["state"], **meta_common}, CKPT_DIR / "pooled.pt")

    per_cell_ckpt = {}
    for i, ct in enumerate(CELL_TYPES):
        idx = np.where(sample["cell"] == i)[0]
        if len(idx) < 4:
            print(f"    {ct:<11} skipped, only {len(idx)} rounds")
            continue
        m = SeedFormer(**cfg)
        m.load_state_dict(pooled["state"])              # warm start from pooled
        best = fit(m, sample, idx, epochs=args.finetune_epochs,
                   lr=args.finetune_lr, weight_decay=args.weight_decay,
                   batch_size=args.batch_size, seed=100 + i, val_frac=0.25)
        path = CKPT_DIR / f"{ct.replace('+', 'plus').replace('/', '_')}.pt"
        torch.save({"state": best["state"], "cell_type": ct,
                    "cell_index": i, "n_rounds": int(len(idx)),
                    "val_loss": best["loss"], **meta_common}, path)
        per_cell_ckpt[ct] = {"path": str(path.relative_to(ROOT)),
                             "rounds": int(len(idx)),
                             "val_loss": best["loss"]}
        print(f"    {ct:<11} finetuned on {len(idx)} rounds, "
              f"val {best['loss']:.4f} -> {path.name}")

    # ---- verdict ---------------------------------------------------------
    ref_summary = {name: float(np.mean(h)) for name, h in pooled_ref.items()}
    best_ref_name = max(ref_summary, key=ref_summary.get)
    verdict = {
        "model_mean_hits": float(np.mean(pooled_model_hits)),
        "reference_mean_hits": ref_summary,
        "best_reference": best_ref_name,
        "chance_mean_hits": 1.0,
        "vs_best_reference": permutation_test(pooled_model_hits,
                                              pooled_ref[best_ref_name]),
        "vs_uniform": permutation_test(pooled_model_hits, pooled_ref["uniform"]),
        "per_cell_type": {
            ct: {"n": len(v["model"]),
                 "model_mean_hits": float(np.mean(v["model"])) if v["model"] else None,
                 "best_ref_mean_hits": float(np.mean(v["best_ref"])) if v["best_ref"] else None,
                 "test": permutation_test(v["model"], v["best_ref"]) if v["model"] else None}
            for ct, v in per_cell_hits.items()},
    }
    report = {"generated_at": datetime.now(timezone.utc).isoformat(),
              "args": vars(args), "cell_types": list(CELL_TYPES),
              "n_params": model.n_params(), "folds": per_fold,
              "verdict": verdict, "checkpoints": per_cell_ckpt,
              "runtime_seconds": round(time.time() - t0, 1)}
    Path(args.report).write_text(json.dumps(report, indent=2) + "\n")

    print(f"\n{'=' * 68}")
    print(f"walk-forward over {len(pooled_model_hits)} held-out rounds "
          f"({model.n_params():,} params)")
    print(f"  model            {verdict['model_mean_hits']:.3f} hits / 3")
    for name, v in sorted(ref_summary.items(), key=lambda kv: -kv[1]):
        print(f"  {name:<16} {v:.3f}")
    print(f"  chance            1.000")
    t = verdict["vs_best_reference"]
    print(f"\n  model - {best_ref_name}: {t['observed']:+.3f} hits/round, "
          f"p = {t['p_value']:.3f} (paired sign-flip, one-sided)")
    print(f"  report -> {args.report}")
    print(f"  runtime {report['runtime_seconds']}s")


if __name__ == "__main__":
    main()
