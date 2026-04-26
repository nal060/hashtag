"""Tests for stamp.encoder: dictionary, field encoding, sequence hash, nonce search."""

from __future__ import annotations

import pytest

from stamp import encoder, hsm_mock


# ---------- dictionary ----------

def test_standard_dictionary_round_trip():
    d = encoder.standard_dictionary()
    bits = [0, 0, 0, 1, 1, 0, 1, 1]  # AC GT
    dna = encoder.bits_to_dna(bits, d)
    assert dna == "ACGT"
    assert encoder.dna_to_bits(dna, d) == bits


def test_generated_dictionary_is_a_permutation():
    seed = encoder.derive_seed(1, 1, 0)
    d = encoder.generate_dictionary(seed)
    assert sorted(d.keys()) == [0, 1, 2, 3]
    assert sorted(d.values()) == ["A", "C", "G", "T"]


def test_generated_dictionary_is_deterministic():
    s1 = encoder.derive_seed(7, 2, 5)
    s2 = encoder.derive_seed(7, 2, 5)
    assert encoder.generate_dictionary(s1) == encoder.generate_dictionary(s2)


def test_dictionary_changes_with_seed():
    seen = set()
    for inc in range(20):
        seed = encoder.derive_seed(1, 1, inc)
        d = encoder.generate_dictionary(seed)
        seen.add(tuple(d[i] for i in range(4)))
    assert len(seen) > 1  # we see multiple permutations


# ---------- bits <-> DNA ----------

def test_bits_to_dna_with_permuted_dictionary():
    seed = encoder.derive_seed(0, 0, 0)
    d = encoder.generate_dictionary(seed)
    bits = [0, 0, 0, 1, 1, 0, 1, 1]
    dna = encoder.bits_to_dna(bits, d)
    assert encoder.dna_to_bits(dna, d) == bits


def test_bits_to_dna_rejects_odd_length():
    with pytest.raises(ValueError):
        encoder.bits_to_dna([0, 1, 0], encoder.standard_dictionary())


# ---------- molecular constraints ----------

def test_constraints_reject_forbidden_site():
    # GCGGCCGC (NotI) is one of many CommOnly sites in the blacklist; its
    # presence guarantees rejection. The exact reported substring depends on
    # blacklist iteration order, so we only assert the rejection class.
    dna = "ACGTACGT" + "GCGGCCGC" + "ACGTACGT" * 3
    ok, reason = encoder.check_molecular_constraints(dna)
    assert not ok
    assert "forbidden site" in reason


def test_constraints_reject_homopolymer():
    # Hand-picked filler with balanced GC and no 6+bp CommOnly site, ending
    # in 5x A so the first failure surfaced is the homopolymer.
    ok, reason = encoder.check_molecular_constraints("GCGAAGTAGTGCAAAAA")
    assert not ok
    assert "homopolymer" in reason


def test_constraints_reject_low_gc():
    dna = "ATATATATAT" * 10
    ok, reason = encoder.check_molecular_constraints(dna)
    assert not ok
    assert "GC" in reason


def test_constraints_accept_clean_dna():
    """Generated DNA that should pass all constraints. We retry until we hit one."""
    import random as _r
    rng = _r.Random(0)
    bases = "ACGT"
    for _ in range(200):
        dna = "".join(rng.choice(bases) for _ in range(50))
        ok, _reason = encoder.check_molecular_constraints(dna)
        if ok:
            return
    pytest.fail("could not generate a clean 50bp dna in 200 tries")


# ---------- sequence hash ----------

def test_sequence_hash_deterministic():
    seq = "ACGT" * 100
    seed = encoder.derive_seed(1, 1, 0)
    a = encoder.compute_sequence_hash(seq, seed)
    b = encoder.compute_sequence_hash(seq, seed)
    assert a == b
    assert len(a) == 2


def test_sequence_hash_changes_with_sequence():
    seed = encoder.derive_seed(1, 1, 0)
    a = encoder.compute_sequence_hash("ACGT" * 100, seed)
    b = encoder.compute_sequence_hash("TGCA" * 100, seed)
    assert a != b


def test_sequence_hash_short_sequence_is_handled():
    seed = encoder.derive_seed(1, 1, 0)
    h = encoder.compute_sequence_hash("AC", seed)
    assert len(h) == 2


# ---------- field block round-trips ----------

def _make_fields(**overrides):
    base = dict(
        increment=0,
        synthesizer_id=42,
        run_counter=7,
        sequence_length=4700,
        mech_hash=b"\xab\xcd",
        chain_hash=b"\x12\x34",
        h_sig=b"\xfe\xed",
        seq_hash=b"\xde\xad",
        landmarks=[0, 1, 2, 3, 0, 1, 2, 3, 0, 1],
    )
    base.update(overrides)
    return encoder.Fields(**base)


def test_plaintext_block_round_trip():
    fields = _make_fields()
    dna = encoder.encode_plaintext_block(fields)
    decoded, syndromes = encoder.decode_plaintext_block(dna)
    assert decoded["increment"] == 0
    assert decoded["synthesizer_id"] == 42
    assert decoded["run_counter"] == 7
    assert decoded["sequence_length"] == 4700
    assert all(s == 0 for s in syndromes)


def test_plaintext_block_corrects_one_bit_error_per_block():
    fields = _make_fields(synthesizer_id=0xABC, run_counter=0x1234)
    dna = encoder.encode_plaintext_block(fields)
    # flip one base in the middle (changes 2 bits — but Hamming corrects 1 per block)
    # so flip just a single bit by changing the corresponding base via dictionary swap
    # easier path: corrupt the 2nd-block DNA at position 0 — only single-bit if mapping aligns
    # we'll simulate a single-bit error by flipping one base to its 1-hamming-neighbor
    bases = list(dna)
    # standard mapping is A=00,C=01,G=10,T=11. A<->C is a single-bit flip in the LSB.
    if bases[0] == "A":
        bases[0] = "C"
    elif bases[0] == "C":
        bases[0] = "A"
    elif bases[0] == "G":
        bases[0] = "T"
    else:
        bases[0] = "G"
    corrupted = "".join(bases)
    decoded, syndromes = encoder.decode_plaintext_block(corrupted)
    assert decoded["synthesizer_id"] == 0xABC
    assert decoded["run_counter"] == 0x1234
    assert any(s != 0 for s in syndromes)


def test_dict_protected_block_round_trip():
    fields = _make_fields()
    seed = encoder.derive_seed(fields.synthesizer_id, fields.run_counter, fields.increment)
    d = encoder.generate_dictionary(seed)
    dna = encoder.encode_dict_protected_block(fields, d)
    decoded, syndromes = encoder.decode_dict_protected_block(dna, d)
    assert decoded["mech_hash"] == fields.mech_hash
    assert decoded["chain_hash"] == fields.chain_hash
    assert decoded["h_sig"] == fields.h_sig
    assert all(s == 0 for s in syndromes)


def test_seq_hash_round_trip():
    fields = _make_fields()
    seed = encoder.derive_seed(1, 1, 0)
    d = encoder.generate_dictionary(seed)
    dna = encoder.encode_seq_hash(fields, d)
    assert encoder.decode_seq_hash(dna, d) == fields.seq_hash


def test_landmarks_round_trip():
    fields = _make_fields(landmarks=[3, 2, 1, 0, 3, 2, 1, 0, 3, 2])
    seed = encoder.derive_seed(1, 1, 0)
    d = encoder.generate_dictionary(seed)
    dna = encoder.encode_landmarks(fields.landmarks, d, seed)
    assert len(dna) == encoder.N_LANDMARKS
    assert encoder.decode_landmarks(dna, d, seed) == fields.landmarks


def test_landmarks_repeating_values_dont_homopolymer_after_seed_shift():
    """The whole point of the seed-derived shift: repeating landmark values
    must not produce a homopolymer that no increment can break.
    """
    repeating = [1] * encoder.N_LANDMARKS
    found_clean = False
    for inc in range(256):
        seed = encoder.derive_seed(1, 1, inc)
        d = encoder.generate_dictionary(seed)
        dna = encoder.encode_landmarks(repeating, d, seed)
        ok, _ = encoder.check_molecular_constraints(dna)
        if ok:
            found_clean = True
            break
    assert found_clean, "could not find a non-homopolymer landmark encoding in 256 increments"


# ---------- nonce search ----------

def _ascii_primer(length: int, seed: int) -> str:
    """Cheap deterministic primer that won't trip constraints on its own."""
    import random as _r
    rng = _r.Random(seed)
    bases = "ACGT"
    while True:
        out = "".join(rng.choice(bases) for _ in range(length))
        ok, _ = encoder.check_molecular_constraints(out)
        if ok:
            return out


def _varied_landmarks(seed: int = 0) -> list[int]:
    import random as _r
    rng = _r.Random(seed)
    return [rng.randrange(4) for _ in range(encoder.N_LANDMARKS)]


def test_find_valid_encoding_succeeds_on_simple_input():
    fwd = _ascii_primer(20, 1)
    rev = _ascii_primer(20, 2)
    sequence = "ACGTACGTGGCATAGCATGCATGGTACGTACGT" * 50
    mech_hash = hsm_mock.compute_mech_hash("idle", "v1.0")
    result = encoder.find_valid_encoding(
        sequence=sequence,
        synthesizer_id=42,
        run_counter=7,
        mech_hash=mech_hash,
        landmarks=_varied_landmarks(),
        primer_fwd=fwd,
        primer_rev=rev,
    )
    assert 0 <= result.increment < 256
    decoded, _ = encoder.decode_plaintext_block(result.plaintext_dna)
    assert decoded["increment"] == result.increment
    assert decoded["synthesizer_id"] == 42
    assert decoded["run_counter"] == 7
    assert decoded["sequence_length"] == len(sequence)


def test_find_valid_encoding_uses_increment_in_seed():
    fwd = _ascii_primer(20, 1)
    rev = _ascii_primer(20, 2)
    sequence = "ACGTACGTGGCATAGCATGCATGGTACGTACGT" * 50
    mech_hash = hsm_mock.compute_mech_hash("idle", "v1.0")
    result = encoder.find_valid_encoding(
        sequence=sequence,
        synthesizer_id=42,
        run_counter=7,
        mech_hash=mech_hash,
        landmarks=_varied_landmarks(),
        primer_fwd=fwd,
        primer_rev=rev,
    )
    seed_zero = encoder.derive_seed(42, 7, 0)
    if result.increment != 0:
        assert encoder.generate_dictionary(seed_zero) != result.dictionary


def test_full_pipeline_returns_signature_for_ledger():
    """The encoder returns the full 32-byte signature so the caller can post
    it to the ledger. h_sig in the barcode is the truncation of SHA256 of it.
    """
    fwd = _ascii_primer(20, 11)
    rev = _ascii_primer(20, 12)
    sequence = "ACGTACGTGGCATAGCATGCATGGTACGTACGTGCATGCATGCATCGTAG" * 30
    mech_hash = hsm_mock.compute_mech_hash("idle", "v1.0")
    result = encoder.find_valid_encoding(
        sequence=sequence,
        synthesizer_id=1,
        run_counter=1,
        mech_hash=mech_hash,
        landmarks=_varied_landmarks(seed=99),
        primer_fwd=fwd,
        primer_rev=rev,
    )
    assert len(result.full_signature) == 32
    assert hsm_mock.h_sig(result.full_signature) == result.fields.h_sig
    body = encoder.assemble_body(result)
    full = fwd + body + rev
    ok, reason = encoder.check_molecular_constraints(full)
    assert ok, reason


def test_pipeline_chain_and_seq_hash_internally_consistent():
    fwd = _ascii_primer(20, 33)
    rev = _ascii_primer(20, 34)
    sequence = "GCATGCATGCATGCAT" * 50
    mech_hash = hsm_mock.compute_mech_hash("running", "v0.9")
    result = encoder.find_valid_encoding(
        sequence=sequence,
        synthesizer_id=99,
        run_counter=3,
        mech_hash=mech_hash,
        landmarks=_varied_landmarks(seed=7),
        primer_fwd=fwd,
        primer_rev=rev,
    )
    expected_seq_hash = encoder.compute_sequence_hash(sequence, result.seed)
    assert result.fields.seq_hash == expected_seq_hash
    assert result.fields.mech_hash == mech_hash
    assert result.fields.chain_hash == hsm_mock.compute_chain_hash(mech_hash, expected_seq_hash)
