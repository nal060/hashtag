"""Mock public ledger for STAMP.

Real deployment would use an append-only log (transparency log, blockchain,
notary service). For the demo: an in-memory dict, optionally persisted to
JSON on disk so the GUI can survive a restart.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Optional


class Ledger:
    def __init__(self, path: Optional[Path] = None) -> None:
        self.path = path
        self._entries: dict[str, dict] = {}
        if path is not None and path.exists():
            self._entries = json.loads(path.read_text())

    @staticmethod
    def _key(synthesizer_id: int, run_counter: int) -> str:
        return f"{synthesizer_id}:{run_counter}"

    def post(
        self,
        synthesizer_id: int,
        run_counter: int,
        signature: bytes,
        timestamp: int,
        seq_hash: bytes,
    ) -> None:
        self._entries[self._key(synthesizer_id, run_counter)] = {
            "signature": signature.hex(),
            "timestamp": timestamp,
            "seq_hash": seq_hash.hex(),
        }
        self._flush()

    def lookup(self, synthesizer_id: int, run_counter: int) -> Optional[dict]:
        entry = self._entries.get(self._key(synthesizer_id, run_counter))
        if entry is None:
            return None
        return {
            "signature": bytes.fromhex(entry["signature"]),
            "timestamp": entry["timestamp"],
            "seq_hash": bytes.fromhex(entry["seq_hash"]),
        }

    def _flush(self) -> None:
        if self.path is not None:
            self.path.write_text(json.dumps(self._entries, indent=2))
