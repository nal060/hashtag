"""Mock public ledger for STAMP.

Real deployment would use an append-only log (transparency log, blockchain,
notary service). For the demo: an in-memory dict, optionally persisted to
JSON on disk so the GUI can survive a restart.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Iterable, Optional


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
        landmark_hits: Optional[Iterable] = None,
    ) -> None:
        """Post an attestation. `landmark_hits` is an optional iterable of
        objects with `slot/site/position/upstream/feature` attributes (e.g.
        landmarks.LandmarkHit); they are serialized as plain dicts so the
        ledger stays free of any landmarks dependency.
        """
        entry: dict = {
            "signature": signature.hex(),
            "timestamp": timestamp,
            "seq_hash": seq_hash.hex(),
        }
        if landmark_hits is not None:
            entry["landmark_hits"] = [
                {"slot": h.slot, "site": h.site, "position": h.position,
                 "upstream": h.upstream, "feature": h.feature}
                for h in landmark_hits
            ]
        self._entries[self._key(synthesizer_id, run_counter)] = entry
        self._flush()

    def lookup(self, synthesizer_id: int, run_counter: int) -> Optional[dict]:
        entry = self._entries.get(self._key(synthesizer_id, run_counter))
        if entry is None:
            return None
        out = {
            "signature": bytes.fromhex(entry["signature"]),
            "timestamp": entry["timestamp"],
            "seq_hash": bytes.fromhex(entry["seq_hash"]),
        }
        if "landmark_hits" in entry:
            out["landmark_hits"] = entry["landmark_hits"]
        return out

    def _flush(self) -> None:
        if self.path is not None:
            self.path.write_text(json.dumps(self._entries, indent=2))
