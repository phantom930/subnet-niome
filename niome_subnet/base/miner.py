# The MIT License (MIT)
# Copyright © 2023 Yuma Rao

# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated
# documentation files (the "Software"), to deal in the Software without restriction, including without limitation
# the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software,
# and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all copies or substantial portions of
# the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO
# THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
# THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
# OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
# DEALINGS IN THE SOFTWARE.

import time
import asyncio
import ipaddress
import os
import socket
import threading
import argparse
import logging
import traceback
import urllib.request

import bittensor as bt
import uvicorn

from fastapi import FastAPI, Request, HTTPException

from niome_subnet.base.neuron import BaseNeuron
from niome_subnet.utils import add_miner_args, fetch_metagraph_with_retry

from typing import Union

logger = logging.getLogger(__name__)


class BaseMinerNeuron(BaseNeuron):
    """
    Base class for Bittensor miners. Uses a FastAPI HTTP server instead of bt.Axon.
    """

    neuron_type: str = "MinerNeuron"

    @classmethod
    def add_args(cls, parser: argparse.ArgumentParser):
        super().add_args(parser)
        add_miner_args(cls, parser)
        # Four separate things, because behind NAT they are genuinely different values: where the
        # server binds locally, and what address:port validators are told to dial.
        parser.add_argument(
            "--axon.ip",
            type=str,
            help="Local address the miner HTTP server binds to.",
            default="0.0.0.0",
        )
        parser.add_argument(
            "--axon.port",
            type=int,
            help="Local port the miner HTTP server binds to.",
            default=8091,
        )
        parser.add_argument(
            "--axon.external-ip",
            "--axon.external_ip",
            dest="axon.external_ip",
            type=str,
            help=(
                "Public IP to publish on chain. Defaults to $AXON_EXTERNAL_IP, then autodiscovery "
                "via a public echo service, then a local hostname lookup. Set this explicitly when "
                "the miner is behind NAT, where the local lookup returns an unroutable address."
            ),
            default=None,
        )
        parser.add_argument(
            "--axon.external-port",
            "--axon.external_port",
            dest="axon.external_port",
            type=int,
            help=(
                "Public port to publish on chain. Defaults to $AXON_EXTERNAL_PORT, then "
                "--axon.port. Set this when the port forward maps a different outside port to "
                "this miner (e.g. 9091 on the router -> 8091 here)."
            ),
            default=None,
        )

    def __init__(self, config=None):
        super().__init__(config=config)

        if not self.config.blacklist.force_validator_permit:
            logger.warning(
                "You are allowing non-validators to send requests to your miner. This is a security risk."
            )
        if self.config.blacklist.allow_non_registered:
            logger.warning(
                "You are allowing non-registered entities to send requests to your miner. This is a security risk."
            )

        axon_config = getattr(self.config, "axon", None)
        self.axon_ip = getattr(axon_config, "ip", None) or "0.0.0.0"
        self.axon_port = getattr(axon_config, "port", 8091)
        self.axon_external_ip = (
            getattr(axon_config, "external_ip", None) or os.getenv("AXON_EXTERNAL_IP")
        )
        self.axon_external_port = int(
            getattr(axon_config, "external_port", None)
            or os.getenv("AXON_EXTERNAL_PORT")
            or self.axon_port
        )

        # Build the FastAPI app; subclass attaches routes in forward/blacklist/priority
        self.app = FastAPI()
        self._setup_routes()

        # Instantiate runners
        self.should_exit: bool = False
        self.is_running: bool = False
        self.thread: Union[threading.Thread, None] = None
        self.lock = asyncio.Lock()

    def _setup_routes(self):
        """Attach the /forward route to the FastAPI app."""
        app = self.app
        miner = self

        @app.post("/forward")
        async def handle_forward(request: Request):
            try:
                headers = dict(request.headers)
                body = await request.body()

                # Verify hotkey-signed request
                try:
                    caller = bt.http_auth.verify(
                        headers,
                        body,
                        method="POST",
                        path="/forward",
                        self_hotkey_ss58=miner.wallet.hotkey.ss58_address,
                    )
                except bt.http_auth.AuthError as e:
                    raise HTTPException(status_code=401, detail=str(e))

                # Run blacklist check
                if await miner.blacklist(caller.hotkey_ss58):
                    raise HTTPException(status_code=403, detail="blacklisted")

                return await miner.forward(body, caller.hotkey_ss58)
            except HTTPException:
                raise
            except Exception as e:
                logger.error(f"Error handling /forward: {e}")
                raise HTTPException(status_code=500, detail=str(e))

    # Echo services used to discover the egress address. Plain-text bodies, queried in order.
    EXTERNAL_IP_SERVICES = (
        "https://checkip.amazonaws.com",
        "https://api.ipify.org",
        "https://ifconfig.me/ip",
    )

    @staticmethod
    def _is_routable(ip: str) -> bool:
        """Whether an address is one a validator on the public internet could dial."""
        try:
            return ipaddress.ip_address(ip).is_global
        except ValueError:
            return False

    def _discover_external_ip(self) -> str | None:
        """Ask a public echo service what address our traffic appears to come from.

        This answers "which IP does the outside world see", which is the right value to advertise
        — but only the *address* half. Whether inbound connections on that address reach this
        process is a property of the NAT in between that no echo service can observe, so a
        discovered address is a starting point for the port forward, never proof of one.
        """
        for url in self.EXTERNAL_IP_SERVICES:
            try:
                with urllib.request.urlopen(url, timeout=5) as response:
                    ip = response.read().decode().strip()
                if self._is_routable(ip):
                    logger.info(f"Discovered external IP {ip} via {url}")
                    return ip
                logger.debug(f"{url} returned {ip!r}, which is not globally routable")
            except Exception as e:
                logger.debug(f"External IP lookup via {url} failed: {e}")
        return None

    def _resolve_external_ip(self) -> str:
        """The address validators will dial, in descending order of authority.

        1. the explicit override, which is the only source that can express "the outside world
           reaches me at an address this machine cannot see" — the normal case behind NAT;
        2. a public echo service, correct whenever egress and ingress share an address;
        3. the local hostname lookup, correct only when this host owns its public address.

        A private result is published rather than rejected — the chain accepts it and some
        operators do run validators on the same private network — but it is the single most
        common reason a healthy-looking miner receives nothing, so it is warned about loudly.
        """
        if self.axon_external_ip:
            return self.axon_external_ip

        discovered = self._discover_external_ip()
        if discovered:
            logger.warning(
                f"Publishing autodiscovered {discovered}:{self.axon_external_port}. That is where "
                "this miner's traffic leaves from, which is not the same as where validators can "
                "reach it — confirm TCP is forwarded to this machine, or pass --axon.external-ip "
                "/ --axon.external-port explicitly."
            )
            return discovered

        ip = socket.gethostbyname(socket.gethostname())
        if not self._is_routable(ip):
            logger.warning(
                f"Publishing {ip}, which is not routable from outside this host — validators will "
                "not be able to reach this miner. Pass --axon.external-ip <PUBLIC_IP> (or set "
                "AXON_EXTERNAL_IP) and forward the port to this machine."
            )
        return ip

    def _serve_axon_on_chain(self):
        """Register this miner's IP:port on chain."""
        try:
            ip = self._resolve_external_ip()
            result = self.subtensor.execute(
                bt.ServeAxon(
                    netuid=self.config.netuid,
                    ip=ip,
                    port=self.axon_external_port,
                ),
                self.wallet,
            )
            # execute() reports chain-level rejections in its return value instead of raising, so
            # the except below never sees them. ServingRateLimitExceeded ("Custom error: 12") is
            # the one that bites on restart, and logging success through it hides the fact that
            # the chain is still pointing validators at the previously published endpoint.
            if result is not None and not getattr(result, "success", True):
                logger.error(
                    f"Chain rejected serve_axon for {ip}:{self.axon_external_port}: "
                    f"{getattr(result, 'message', result)}. Validators will keep using whatever "
                    "endpoint was published before this."
                )
                return
            logger.info(
                f"Served miner axon on network: {self.config.network} netuid: "
                f"{self.config.netuid} | bound {self.axon_ip}:{self.axon_port} -> advertised "
                f"{ip}:{self.axon_external_port}"
            )
        except Exception as e:
            logger.error(f"Failed to serve axon on chain: {e}")

    def run(self):
        """
        Initiates and manages the main loop for the miner on the Bittensor network.
        """

        # Check that miner is registered on the network.
        self.sync()

        # Register axon on chain.
        self._serve_axon_on_chain()

        logger.info(f"Miner starting at block: {self.block}")

        # Start the FastAPI server in a daemon thread.
        server_config = uvicorn.Config(
            self.app,
            host=self.axon_ip,
            port=self.axon_port,
            log_level="warning",
        )
        server = uvicorn.Server(server_config)
        server_thread = threading.Thread(target=server.run, daemon=True)
        server_thread.start()

        # This loop maintains the miner's operations until intentionally stopped.
        try:
            while not self.should_exit:
                try:
                    # Wait out the epoch instead of spinning. sync() is two chain calls — a
                    # registration lookup and a full metagraph fetch — and nothing the miner does
                    # needs to react sooner than the epoch_length the operator configured. The
                    # one-second granularity is only so should_exit stays responsive; the block
                    # height behind should_sync_metagraph is itself TTL-cached, so polling it
                    # this often costs nothing.
                    while not self.should_exit and not self.should_sync_metagraph():
                        time.sleep(1)
                    if self.should_exit:
                        break

                    self.sync()
                    self.step += 1

                except Exception as err:
                    logger.warning("Miner step error (will retry): %s", err)
                    logger.debug(traceback.format_exc())
                    time.sleep(12)

        except KeyboardInterrupt:
            logger.info("Miner killed by keyboard interrupt.")
            server.should_exit = True
            exit()
        finally:
            server.should_exit = True

    def run_in_background_thread(self):
        """Starts the miner's operations in a separate background thread."""
        if not self.is_running:
            logger.debug("Starting miner in background thread.")
            self.should_exit = False
            self.thread = threading.Thread(target=self.run, daemon=True)
            self.thread.start()
            self.is_running = True
            logger.debug("Started")

    def stop_run_thread(self):
        """Stops the miner's operations that are running in the background thread."""
        if self.is_running:
            logger.debug("Stopping miner in background thread.")
            self.should_exit = True
            if self.thread is not None:
                self.thread.join(5)
            self.is_running = False
            logger.debug("Stopped")

    def __enter__(self):
        self.run_in_background_thread()
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.stop_run_thread()

    def resync_metagraph(self):
        """Resyncs the metagraph."""
        logger.info("resync_metagraph()")
        self.metagraph = fetch_metagraph_with_retry(self.subtensor, self.netuid, commitments=False)
