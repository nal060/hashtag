"""Tests for hsm_mock and ledger."""

from __future__ import annotations

import json
from pathlib import Path

from stamp import hsm_mock, ledger as ledger_mod


def test_mock_hsm_sign_is_deterministic():
    a = hsm_mock.mock_hsm_sign(b"\xaa\xbb", b"\xab\xcd", b"\x12\x34")
    b = hsm_mock.mock_hsm_sign(b"\xaa\xbb", b"\xab\xcd", b"\x12\x34")
    assert a == b
    assert len(a) == 32


def test_mock_hsm_sign_changes_with_each_hash():
    """Changing any of the three signed hashes must change the signature
    — that is the whole point of binding mech/seq/chain together."""
    base = hsm_mock.mock_hsm_sign(b"\x00\x00", b"\x00\x00", b"\x00\x00")
    diff_mech = hsm_mock.mock_hsm_sign(b"\x01\x00", b"\x00\x00", b"\x00\x00")
    diff_seq = hsm_mock.mock_hsm_sign(b"\x00\x00", b"\x01\x00", b"\x00\x00")
    diff_chain = hsm_mock.mock_hsm_sign(b"\x00\x00", b"\x00\x00", b"\x01\x00")
    assert len({base, diff_mech, diff_seq, diff_chain}) == 4


def test_h_sig_is_16_bit_truncation():
    sig = hsm_mock.mock_hsm_sign(b"\x00\x00", b"\x00\x00", b"\x00\x00")
    h = hsm_mock.h_sig(sig)
    assert len(h) == 2


def test_mech_and_chain_hashes():
    mech = hsm_mock.compute_mech_hash("idle", "v1.2.3")
    seq = b"\xde\xad"
    chain = hsm_mock.compute_chain_hash(mech, seq)
    assert len(mech) == 2
    assert len(chain) == 2
    # changing either input changes the chain
    chain2 = hsm_mock.compute_chain_hash(mech, b"\xbe\xef")
    assert chain != chain2


def test_ledger_post_and_lookup_round_trip():
    L = ledger_mod.Ledger()
    sig = b"\x01" * 32
    seq_hash = b"\x42\x43"
    L.post(synthesizer_id=7, run_counter=3, signature=sig, timestamp=1700000000, seq_hash=seq_hash)
    entry = L.lookup(7, 3)
    assert entry is not None
    assert entry["signature"] == sig
    assert entry["seq_hash"] == seq_hash
    assert entry["timestamp"] == 1700000000


def test_ledger_lookup_missing_returns_none():
    L = ledger_mod.Ledger()
    assert L.lookup(0, 0) is None


def test_ledger_persists_to_disk(tmp_path: Path):
    path = tmp_path / "ledger.json"
    L1 = ledger_mod.Ledger(path)
    L1.post(1, 0, b"\xaa" * 32, 1700000000, b"\x12\x34")
    # new instance reads the same file
    L2 = ledger_mod.Ledger(path)
    entry = L2.lookup(1, 0)
    assert entry is not None
    assert entry["signature"] == b"\xaa" * 32
    # verify the file is well-formed json
    assert "1:0" in json.loads(path.read_text())
