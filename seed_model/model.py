"""SeedFormer: a small pre-norm Transformer over the round history.

Deliberately modest. With 22-32 rounds per cell type the binding constraint is
data, not capacity, so the architecture buys its modernity in regularisation
and optimisation rather than in size:

  * pre-norm blocks with RMSNorm - stable at tiny batch sizes, no bias terms
  * SwiGLU feed-forward - better parameter efficiency than GELU-MLP
  * rotary position embeddings - relative positions, no learned table to overfit
  * masked attention over the left-padded history
  * attention-pooled readout rather than last-token, so a short history is not
    dominated by whichever round happens to sit at the end
  * a cell-type embedding, so one model can pretrain on every cell type and
    then be finetuned per cell type
  * stochastic depth + dropout + weight decay, and an EMA shadow used for eval

The head emits nine logits read as a distribution over the class of a single
draw. The target is the next round's three draws normalised to sum to one, so
a class drawn twice carries twice the mass. That keeps the loss a plain soft
cross-entropy and lets a prediction repeat a class, which real triples do.
"""

import math

import torch
import torch.nn as nn
import torch.nn.functional as F

from .data import CELL_TYPES, FEATURE_DIM, N_CLASSES


class RMSNorm(nn.Module):
    def __init__(self, dim, eps=1e-6):
        super().__init__()
        self.weight = nn.Parameter(torch.ones(dim))
        self.eps = eps

    def forward(self, x):
        return self.weight * x * torch.rsqrt(x.pow(2).mean(-1, keepdim=True) + self.eps)


def rope_cache(seq, dim, device, base=10000.0):
    half = dim // 2
    freqs = torch.exp(-math.log(base) * torch.arange(half, device=device) / half)
    pos = torch.arange(seq, device=device, dtype=torch.float32)
    ang = pos[:, None] * freqs[None, :]
    return torch.cos(ang), torch.sin(ang)


def apply_rope(x, cos, sin):
    """x: (b, heads, seq, head_dim) - rotate consecutive pairs."""
    x1, x2 = x[..., 0::2], x[..., 1::2]
    c, s = cos[None, None], sin[None, None]
    return torch.stack([x1 * c - x2 * s, x1 * s + x2 * c], dim=-1).flatten(-2)


class Attention(nn.Module):
    def __init__(self, dim, heads, dropout):
        super().__init__()
        assert dim % heads == 0
        self.heads, self.head_dim = heads, dim // heads
        self.qkv = nn.Linear(dim, 3 * dim, bias=False)
        self.proj = nn.Linear(dim, dim, bias=False)
        self.dropout = dropout

    def forward(self, x, mask, cos, sin):
        b, t, d = x.shape
        q, k, v = self.qkv(x).chunk(3, dim=-1)
        shape = (b, t, self.heads, self.head_dim)
        q, k, v = (z.view(shape).transpose(1, 2) for z in (q, k, v))
        q, k = apply_rope(q, cos, sin), apply_rope(k, cos, sin)
        # mask is 1 on real rounds; block attention *into* padding
        attn_mask = mask[:, None, None, :].bool().expand(b, 1, t, t)
        out = F.scaled_dot_product_attention(
            q, k, v, attn_mask=attn_mask,
            dropout_p=self.dropout if self.training else 0.0)
        return self.proj(out.transpose(1, 2).reshape(b, t, d))


class SwiGLU(nn.Module):
    def __init__(self, dim, mult=8 / 3):
        super().__init__()
        hidden = int(dim * mult / 2) * 2
        self.w12 = nn.Linear(dim, 2 * hidden, bias=False)
        self.w3 = nn.Linear(hidden, dim, bias=False)

    def forward(self, x):
        a, b = self.w12(x).chunk(2, dim=-1)
        return self.w3(F.silu(a) * b)


class Block(nn.Module):
    def __init__(self, dim, heads, dropout, drop_path):
        super().__init__()
        self.n1, self.n2 = RMSNorm(dim), RMSNorm(dim)
        self.attn = Attention(dim, heads, dropout)
        self.ff = SwiGLU(dim)
        self.drop = nn.Dropout(dropout)
        self.drop_path = drop_path

    def _maybe_drop(self, residual):
        if not self.training or self.drop_path <= 0:
            return residual
        keep = 1.0 - self.drop_path
        gate = torch.rand(residual.shape[0], 1, 1, device=residual.device) < keep
        return residual * gate / keep

    def forward(self, x, mask, cos, sin):
        x = x + self._maybe_drop(self.drop(self.attn(self.n1(x), mask, cos, sin)))
        return x + self._maybe_drop(self.drop(self.ff(self.n2(x))))


class SeedFormer(nn.Module):
    def __init__(self, dim=64, depth=3, heads=4, dropout=0.2, drop_path=0.1,
                 n_cells=len(CELL_TYPES)):
        super().__init__()
        self.dim = dim
        self.input = nn.Sequential(nn.Linear(FEATURE_DIM, dim), RMSNorm(dim))
        self.cell = nn.Embedding(n_cells, dim)
        nn.init.normal_(self.cell.weight, std=0.02)
        self.blocks = nn.ModuleList([
            Block(dim, heads, dropout, drop_path * i / max(depth - 1, 1))
            for i in range(depth)])
        self.norm = RMSNorm(dim)
        self.pool = nn.Linear(dim, 1, bias=False)          # attention pooling
        self.head = nn.Linear(dim, N_CLASSES)
        nn.init.zeros_(self.head.bias)
        nn.init.normal_(self.head.weight, std=0.01)        # start near-uniform

    def forward(self, x, mask, cell):
        h = self.input(x) + self.cell(cell)[:, None, :]
        cos, sin = rope_cache(x.shape[1], self.dim // self.blocks[0].attn.heads,
                              x.device)
        for blk in self.blocks:
            h = blk(h, mask, cos, sin)
        h = self.norm(h)
        w = self.pool(h).squeeze(-1).masked_fill(mask == 0, float("-inf"))
        h = (h * torch.softmax(w, dim=-1)[..., None]).sum(1)
        return self.head(h)

    def n_params(self):
        return sum(p.numel() for p in self.parameters())


class EMA:
    """Shadow weights; evaluation uses these, which is worth real accuracy
    when the training set is a few dozen samples."""

    def __init__(self, model, decay=0.995):
        self.decay = decay
        self.shadow = {k: v.detach().clone() for k, v in model.state_dict().items()}

    @torch.no_grad()
    def update(self, model):
        for k, v in model.state_dict().items():
            if v.dtype.is_floating_point:
                self.shadow[k].mul_(self.decay).add_(v.detach(), alpha=1 - self.decay)
            else:
                self.shadow[k].copy_(v)

    def copy_to(self, model):
        model.load_state_dict(self.shadow)


def soft_cross_entropy(logits, target, smoothing=0.05):
    """CE against a soft target multiset, with uniform label smoothing."""
    if smoothing > 0:
        target = (1 - smoothing) * target + smoothing / target.shape[-1]
    return -(target * F.log_softmax(logits, dim=-1)).sum(-1).mean()
