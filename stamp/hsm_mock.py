"""Mock HSM signing and hash composition for STAMP.

Real HSMs would expose a sign() that wraps a private key. For the demo we use
a deterministic PRG seeded by the inputs so test cases are reproducible and
no key material has to live in the repo.
"""

from __future__ import annotations

import hashlib
import random


# Hash truncation length, in bytes. Must stay in sync with encoder.py's
# H_SIG_BITS / MECH_HASH_BITS / CHAIN_HASH_BITS (all 8 * _HASH_BYTES).
_HASH_BYTES = 3


def mock_hsm_sign(
    mech_hash: bytes,
    seq_hash: bytes,
    chain_hash: bytes,
    seed: int = 42,
) -> bytes:
    """Return a 32-byte deterministic 'signature' over (mech, seq, chain).

    The signature payload is `mech_hash || seq_hash || chain_hash`, binding all
    three barcode hashes together so that a consistent triple-replacement at
    the molecular level still breaks H_sig. (In real deployment this would be
    a real PKI signature under the synthesizer's HSM private key; the mock
    here is a deterministic public function so the demo is reproducible —
    that's a known mock limitation, not a design flaw.)
    """
    rng = random.Random(
        f"{seed}|{mech_hash.hex()}|{seq_hash.hex()}|{chain_hash.hex()}"
    )
    return bytes(rng.randint(0, 255) for _ in range(32))


def h_sig(signature: bytes) -> bytes:
    """24-bit truncation of SHA-256(signature) — the value stored in the barcode."""
    return hashlib.sha256(signature).digest()[:_HASH_BYTES]


def compute_mech_hash(machine_state: str, firmware_version: str) -> bytes:
    """24-bit truncation of SHA-256(machine_state || firmware_version)."""
    payload = f"{machine_state}||{firmware_version}".encode()
    return hashlib.sha256(payload).digest()[:_HASH_BYTES]


def compute_chain_hash(mech_hash: bytes, seq_hash: bytes) -> bytes:
    """24-bit truncation of SHA-256(mech_hash || seq_hash)."""
    return hashlib.sha256(mech_hash + seq_hash).digest()[:_HASH_BYTES]
