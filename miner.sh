#!/usr/bin/env bash
set -euo pipefail

# Run from the project root so `-m` puts it on sys.path and niome_subnet resolves.
cd "$(dirname "${BASH_SOURCE[0]}")"

python -m neurons.miner --netuid 55 --wallet phantom_coldkey --wallet-hotkey hotkey-1 --axon.external_ip 184.144.255.144 --axon.external_port 52384 --axon.ip 0.0.0.0 --axon.port 9002
