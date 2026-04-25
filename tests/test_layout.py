"""Tests for stamp.layout: barcode assembly + primer-stripping extraction."""

from __future__ import annotations

import random as _r

import pytest

from stamp import encoder, hsm_mock, layout as layout_mod


def _primer(length: int, seed: int) -> str:
    rng = _r.Random(seed)
    while True:
        out = "".join(rng.choice("ACGT") for _ in range(length))
        if encoder.check_molecular_constraints(out)[0]:
            return out


def _varied_landmarks(seed: int = 0) -> list[int]:
    rng = _r.Random(seed)
    return [rng.randrange(4) for _ in range(encoder.N_LANDMARKS)]


def _make_result(layout: layout_mod.BarcodeLayout, *, sequence: str | None = None):
    if sequence is None:
        sequence = "ACGTACGTGGCATAGCATGCATGGTACGTACGTGCATGCATGCATCGTAG" * 30
    mech = hsm_mock.compute_mech_hash("idle", "v1.0")
    return encoder.find_valid_encoding(
        sequence=sequence,
        synthesizer_id=42,
        run_counter=7,
        mech_hash=mech,
        landmarks=_varied_landmarks(seed=5),
        primer_fwd=layout.primer_fwd,
        primer_rev=layout.primer_rev,
        layout=layout,
    )


def test_layout_lengths_are_consistent():
    fwd = _primer(20, 100); rev = _primer(20, 101)
    layout = layout_mod.BarcodeLayout(primer_fwd=fwd, primer_rev=rev)
    assert layout.plaintext_len == 31
    assert layout.dict_protected_len == 31
    assert layout.seq_hash_len == 8
    assert layout.body_len == 31 + 31 + 8 + 10  # plaintext + dict + seqhash + 10 landmarks
    assert layout.full_len == layout.body_len + 40


def test_layout_invalid_cluster_sizes_raise():
    fwd = _primer(20, 100); rev = _primer(20, 101)
    with pytest.raises(ValueError):
        layout_mod.BarcodeLayout(
            primer_fwd=fwd, primer_rev=rev,
            landmarks_after_fwd_primer=2,
            landmarks_between_plaintext_and_dict=1,
            landmarks_between_dict_and_seqhash=1,
            landmarks_before_rev_primer=5,  # 2+1+1+5 = 9, not 10
        )


def test_assemble_produces_correct_total_length():
    fwd = _primer(20, 100); rev = _primer(20, 101)
    layout = layout_mod.BarcodeLayout(primer_fwd=fwd, primer_rev=rev)
    result = _make_result(layout)
    body = layout.assemble(result)
    assert len(body) == layout.body_len
    full = layout.assemble_full(result)
    assert len(full) == layout.full_len
    assert full.startswith(fwd)
    assert full.endswith(rev)


def test_extract_recovers_original_blocks():
    fwd = _primer(20, 200); rev = _primer(20, 201)
    layout = layout_mod.BarcodeLayout(primer_fwd=fwd, primer_rev=rev)
    result = _make_result(layout)
    body = layout.assemble(result)
    parts = layout.extract(body)
    assert parts["plaintext_dna"] == result.plaintext_dna
    assert parts["dict_protected_dna"] == result.dict_protected_dna
    assert parts["seq_hash_dna"] == result.seq_hash_dna
    assert parts["landmarks_dna"] == result.landmarks_dna


def test_extract_full_strips_primers():
    fwd = _primer(20, 200); rev = _primer(20, 201)
    layout = layout_mod.BarcodeLayout(primer_fwd=fwd, primer_rev=rev)
    result = _make_result(layout)
    full = layout.assemble_full(result)
    parts = layout.extract_full(full)
    assert parts["plaintext_dna"] == result.plaintext_dna


def test_extract_full_rejects_primer_mismatch():
    fwd = _primer(20, 200); rev = _primer(20, 201)
    layout = layout_mod.BarcodeLayout(primer_fwd=fwd, primer_rev=rev)
    result = _make_result(layout)
    full = layout.assemble_full(result)
    bad = "AAAAAAAAAAAAAAAAAAAA" + full[20:]
    with pytest.raises(ValueError, match="forward primer"):
        layout.extract_full(bad)


def test_assembled_barcode_passes_constraints():
    """The encoder's increment search uses the layout, so the assembled
    barcode must clear molecular constraints by construction."""
    fwd = _primer(20, 300); rev = _primer(20, 301)
    layout = layout_mod.BarcodeLayout(primer_fwd=fwd, primer_rev=rev)
    result = _make_result(layout)
    full = layout.assemble_full(result)
    ok, reason = encoder.check_molecular_constraints(full)
    assert ok, reason
