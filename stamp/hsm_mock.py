"""Mock HSM signing and hash composition for STAMP.

Real HSMs would expose a sign() that wraps a private key. For the demo we use
a deterministic PRG seeded by the inputs so test cases are reproducible and
no key material has to live in the repo.
"""

from __future__ import annotations

import hashlib
import random


def mock_hsm_sign(
    seq_hash: bytes,
    synthesizer_id: int,
    run_counter: int,
    seed: int = 42,
) -> bytes:
    """Return a 32-byte deterministic 'signature' over the inputs.

    Spec note: signature = mock_hsm_sign(seq_hash || synthesizer_id || run_counter).
    Fixed `seed` keeps test cases reproducible across runs.
    """
    rng = random.Random(f"{seed}{seq_hash.hex()}{synthesizer_id}{run_counter}")
    return bytes(rng.randint(0, 255) for _ in range(32))


def h_sig(signature: bytes) -> bytes:
    """16-bit truncation of SHA-256(signature) — the value stored in the barcode."""
    return hashlib.sha256(signature).digest()[:2]


def compute_mech_hash(machine_state: str, firmware_version: str) -> bytes:
    """16-bit truncation of SHA-256(machine_state || firmware_version)."""
    payload = f"{machine_state}||{firmware_version}".encode()
    return hashlib.sha256(payload).digest()[:2]


def compute_chain_hash(mech_hash: bytes, seq_hash: bytes) -> bytes:
    """16-bit truncation of SHA-256(mech_hash || seq_hash)."""
    return hashlib.sha256(mech_hash + seq_hash).digest()[:2]
