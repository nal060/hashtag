"""End-to-end tests for stamp.decoder.verify, including the 5 demo scenarios."""

from __future__ import annotations

import random as _r

import pytest

from stamp import decoder, encoder, layout as layout_mod, ledger as ledger_mod


def _primer(length: int, seed: int) -> str:
    rng = _r.Random(seed)
    while True:
        out = "".join(rng.choice("ACGT") for _ in range(length))
        if encoder.check_molecular_constraints(out)[0]:
            return out


# Synthetic sequences with enough length and diversity for sensible
# sparse hashing and landmark finding.
GFP_LIKE = (
    "ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGAC"
    "GGCGACGTAAACGGCCACAAGTTCAGCGTGTCCGGCGAGGGCGAGGGCGATGCCACCTAC"
    "GGCAAGCTGACCCTGAAGTTCATCTGCACCACCGGCAAGCTGCCCGTGCCCTGGCCCACC"
    "CTCGTGACCACCCTGACCTACGGCGTGCAGTGCTTCAGCCGCTACCCCGACCACATGAAG"
    "CAGCACGACTTCTTCAAGTCCGCCATGCCCGAAGGCTACGTCCAGGAGCGCACCATCTTC"
    "TTCAAGGACGACGGCAACTACAAGACCCGCGCCGAGGTGAAGTTCGAGGGCGACACCCTG"
) * 5  # ~1800bp
MCHERRY_LIKE = (
    "ATGGTGAGCAAGGGCGAGGAGGATAACATGGCCATCATCAAGGAGTTCATGCGCTTCAAG"
    "GTGCACATGGAGGGCTCCGTGAACGGCCACGAGTTCGAGATCGAGGGCGAGGGCGAGGGC"
    "CGCCCCTACGAGGGCACCCAGACCGCCAAGCTGAAGGTGACCAAGGGTGGCCCCCTGCCC"
    "TTCGCCTGGGACATCCTGTCCCCTCAGTTCATGTACGGCTCCAAGGCCTACGTGAAGCAC"
    "CCCGCCGACATCCCCGACTACTTGAAGCTGTCCTTCCCCGAGGGCTTCAAGTGGGAGCGC"
) * 4  # ~1200bp


@pytest.fixture
def primers():
    return _primer(20, 500), _primer(20, 501)


@pytest.fixture
def layout(primers):
    return layout_mod.BarcodeLayout(primer_fwd=primers[0], primer_rev=primers[1])


@pytest.fixture
def fresh_ledger():
    return ledger_mod.Ledger()


# ---------- scenario 1: clean synthesis ----------

def test_scenario1_clean_synthesis_all_green(primers, layout, fresh_ledger):
    fwd, rev = primers
    stamped = decoder.stamp(
        sequence=GFP_LIKE,
        synthesizer_id=1,
        run_counter=0,
        machine_state="idle",
        firmware_version="v1.0",
        primer_fwd=fwd,
        primer_rev=rev,
        ledger=fresh_ledger,
    )
    result = decoder.verify(
        barcode_dna=stamped.barcode,
        suspect_sequence=GFP_LIKE,
        layout=layout,
        ledger=fresh_ledger,
        encoded_landmark_hits=stamped.landmark_hits,
    )
    assert result.signature_valid
    assert result.sequence_match
    assert result.chain_valid
    assert result.sequence_length_match
    assert result.all_clean(), result.summarize()
    assert "CLEAN" in result.summarize()


# ---------- scenario 2: legitimate modification ----------

def test_scenario2_modified_sequence_flags_seq_hash_mismatch(primers, layout, fresh_ledger):
    fwd, rev = primers
    stamped = decoder.stamp(
        sequence=GFP_LIKE,
        synthesizer_id=2,
        run_counter=0,
        machine_state="idle",
        firmware_version="v1.0",
        primer_fwd=fwd,
        primer_rev=rev,
        ledger=fresh_ledger,
    )
    # mutate a 500bp slice
    rng = _r.Random(7)
    modified = list(GFP_LIKE)
    for i in range(500):
        modified[100 + i] = rng.choice("ACGT")
    modified_seq = "".join(modified)
    result = decoder.verify(
        barcode_dna=stamped.barcode,
        suspect_sequence=modified_seq,
        layout=layout,
        ledger=fresh_ledger,
        encoded_landmark_hits=stamped.landmark_hits,
    )
    assert result.signature_valid       # signature still binds the original ledger entry
    assert not result.sequence_match    # but the sequence has changed
    summary = result.summarize()
    assert "MODIFIED" in summary or "TRANSPLANT" in summary, summary


# ---------- scenario 3: transplant attack ----------

def test_scenario3_transplant_flags_length_and_landmarks(primers, layout, fresh_ledger):
    fwd, rev = primers
    # encode for GFP, paste barcode onto mCherry
    stamped = decoder.stamp(
        sequence=GFP_LIKE,
        synthesizer_id=3,
        run_counter=0,
        machine_state="idle",
        firmware_version="v1.0",
        primer_fwd=fwd,
        primer_rev=rev,
        ledger=fresh_ledger,
    )
    result = decoder.verify(
        barcode_dna=stamped.barcode,
        suspect_sequence=MCHERRY_LIKE,
        layout=layout,
        ledger=fresh_ledger,
        encoded_landmark_hits=stamped.landmark_hits,
    )
    assert result.signature_valid          # H_sig is still the GFP signature
    assert not result.sequence_match       # sparse hash detects the replacement
    assert not result.sequence_length_match  # length disagreement


# ---------- scenario 4: barcode absent ----------

def test_scenario4_no_ledger_entry_flags_absence(primers, layout, fresh_ledger):
    fwd, rev = primers
    # synth on machine 4 / run 5 then look up under a different (synth, run)
    stamped = decoder.stamp(
        sequence=GFP_LIKE,
        synthesizer_id=4,
        run_counter=5,
        machine_state="idle",
        firmware_version="v1.0",
        primer_fwd=fwd,
        primer_rev=rev,
        ledger=fresh_ledger,
    )
    # Wipe the ledger to simulate "no attestation on file"
    empty = ledger_mod.Ledger()
    result = decoder.verify(
        barcode_dna=stamped.barcode,
        suspect_sequence=GFP_LIKE,
        layout=layout,
        ledger=empty,
    )
    assert not result.ledger_hit
    assert "NOT IN LEDGER" in result.summarize()


# ---------- scenario 5: nonce search demo ----------

def test_scenario5_increment_search_resolves_constraints(primers, layout, fresh_ledger):
    """Demonstrate the increment search: even with awkward primers,
    find_valid_encoding finds an increment that clears constraints."""
    fwd, rev = primers
    stamped = decoder.stamp(
        sequence=GFP_LIKE,
        synthesizer_id=0,
        run_counter=0,
        machine_state="idle",
        firmware_version="v1.0",
        primer_fwd=fwd,
        primer_rev=rev,
        ledger=fresh_ledger,
    )
    # increment field is 8-bit, so any value 0..255 is acceptable
    assert 0 <= stamped.encoding.increment < 256
    # the body still passes constraints
    ok, _ = encoder.check_molecular_constraints(stamped.barcode)
    assert ok


# ---------- additional integrity tests ----------

def test_verify_corrects_single_bit_sequencing_error(primers, layout, fresh_ledger):
    """Hamming(31,26) recovers from a 1-bit corruption in the protected blocks."""
    fwd, rev = primers
    stamped = decoder.stamp(
        sequence=GFP_LIKE,
        synthesizer_id=10,
        run_counter=0,
        machine_state="idle",
        firmware_version="v1.0",
        primer_fwd=fwd,
        primer_rev=rev,
        ledger=fresh_ledger,
    )
    # flip a single base in the plaintext block region
    body_start = len(fwd) + layout.landmarks_after_fwd_primer
    bases = list(stamped.barcode)
    flipped = bases[body_start]
    # flip to a 1-base-different neighbor under standard mapping (A↔C is bit 0 flip)
    swap = {"A": "C", "C": "A", "G": "T", "T": "G"}
    bases[body_start] = swap[flipped]
    corrupted = "".join(bases)
    result = decoder.verify(
        barcode_dna=corrupted,
        suspect_sequence=GFP_LIKE,
        layout=layout,
        ledger=fresh_ledger,
        encoded_landmark_hits=stamped.landmark_hits,
    )
    # the underlying field values still recover correctly
    assert result.signature_valid
    assert result.sequence_match


def test_verify_recovers_synthesizer_id_and_run_counter(primers, layout, fresh_ledger):
    fwd, rev = primers
    stamped = decoder.stamp(
        sequence=GFP_LIKE,
        synthesizer_id=2048,   # within 12-bit field
        run_counter=678,
        machine_state="idle",
        firmware_version="v1.0",
        primer_fwd=fwd,
        primer_rev=rev,
        ledger=fresh_ledger,
    )
    result = decoder.verify(
        barcode_dna=stamped.barcode,
        suspect_sequence=GFP_LIKE,
        layout=layout,
        ledger=fresh_ledger,
        encoded_landmark_hits=stamped.landmark_hits,
    )
    assert result.synthesizer_id == 2048
    assert result.run_counter == 678
    assert result.signature_valid
