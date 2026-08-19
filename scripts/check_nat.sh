#!/usr/bin/env bash
# check_nat.sh
#
# Answers one question: does a TCP connection opened from the public internet actually reach
# this machine's miner port?
#
# Testing that from here is not as simple as curl'ing our own public address — most NATs do not
# hairpin, so a local attempt fails even when inbound forwarding is perfectly configured. The
# only trustworthy method is to have a third party connect *inward* while we watch the wire, so
# this drives public fetch services at our advertised endpoint and captures what arrives.
#
# Usage:
#   sudo ./scripts/check_nat.sh [PORT] [EXTERNAL_IP]
#
# Requires tcpdump (packet capture must run as root).

set -uo pipefail

PORT="${1:-8091}"
EXTERNAL_IP="${2:-}"

if [[ $EUID -ne 0 ]]; then
    echo "ERROR: needs root for tcpdump. Re-run with sudo." >&2
    exit 1
fi
if ! command -v tcpdump &>/dev/null; then
    echo "ERROR: tcpdump not installed (apt-get install -y tcpdump)." >&2
    exit 1
fi

if [[ -z "$EXTERNAL_IP" ]]; then
    for svc in https://checkip.amazonaws.com https://api.ipify.org https://ifconfig.me/ip; do
        EXTERNAL_IP="$(curl -s --max-time 5 "$svc" | tr -d '[:space:]')"
        [[ -n "$EXTERNAL_IP" ]] && break
    done
fi
[[ -z "$EXTERNAL_IP" ]] && { echo "ERROR: could not determine external IP; pass it explicitly." >&2; exit 1; }

echo "Target: ${EXTERNAL_IP}:${PORT}"
echo

# ---------------------------------------------------------------------------
# 1. Local preconditions. A failure here is not a NAT problem.
# ---------------------------------------------------------------------------

echo "[1/3] local listener"
LISTEN="$(ss -tlnp 2>/dev/null | grep ":${PORT} ")"
if [[ -z "$LISTEN" ]]; then
    echo "  FAIL nothing is listening on ${PORT} — start the miner first."
    exit 1
fi
echo "  ok  $(echo "$LISTEN" | awk '{print $4}')"
case "$LISTEN" in
    *"0.0.0.0:${PORT}"*|*"*:${PORT}"*|*"[::]:${PORT}"*) ;;
    *) echo "  WARN bound to a specific address, not 0.0.0.0 — forwarded packets may be dropped." ;;
esac

echo "[2/3] host firewall"
if command -v ufw &>/dev/null && ufw status 2>/dev/null | grep -q "^Status: active"; then
    if ufw status | grep -q "${PORT}"; then
        echo "  ok  ufw active and allows ${PORT}"
    else
        echo "  FAIL ufw is active and has no rule for ${PORT}: sudo ufw allow ${PORT}/tcp"
    fi
else
    echo "  ok  ufw inactive or absent"
fi

# ---------------------------------------------------------------------------
# 2. The actual test: watch the wire while third parties connect inward.
# ---------------------------------------------------------------------------

echo "[3/3] inbound reachability"
CAP="$(mktemp)"
tcpdump -i any -n "tcp port ${PORT}" > "$CAP" 2>/dev/null &
TCPDUMP_PID=$!
trap 'kill $TCPDUMP_PID 2>/dev/null; rm -f "$CAP"' EXIT
sleep 3

TARGET="http://${EXTERNAL_IP}:${PORT}/forward"
for prober in \
    "https://api.codetabs.com/v1/proxy/?quest=${TARGET}" \
    "https://api.allorigins.win/raw?url=${TARGET}" \
    "https://r.jina.ai/${TARGET}" ; do
    printf '  probing via %-20s ' "$(echo "$prober" | cut -d/ -f3)"
    timeout 25 curl -s -o /dev/null -w "(prober said HTTP %{http_code})\n" "$prober" || echo "(prober timed out)"
done

sleep 4
kill $TCPDUMP_PID 2>/dev/null
sleep 1

# Only packets from off-box count. Loopback and the local subnet prove nothing about NAT.
# grep -c prints its count and exits 1 when that count is zero, so the failure has to be
# swallowed without appending a second number to the output.
INBOUND="$(grep -c "Flags \[S\]" "$CAP" 2>/dev/null || true)"
INBOUND="${INBOUND:-0}"
echo
if [[ "$INBOUND" -gt 0 ]]; then
    echo "PASS  ${INBOUND} inbound SYN(s) reached this host — forwarding works."
    grep "Flags \[S\]" "$CAP" | head -3
else
    echo "FAIL  no inbound packet arrived. The prober's connection never got here, so the"
    echo "      break is upstream of this machine, not in the miner or its firewall:"
    echo "        - hypervisor: QEMU user-mode needs hostfwd=tcp::${PORT}-:${PORT}, or switch to bridged networking"
    echo "        - router:     forward TCP ${PORT} to this host's LAN address"
    echo "        - ISP:        if the router's WAN address is in 100.64.0.0/10 you are behind CGNAT"
    echo "                      and no forward can work — use a relay host with a public IP"
fi
