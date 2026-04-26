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


def test_verify_detects_plaintext_length_tamper_via_ledger(primers, layout, fresh_ledger):
    """Attacker rewrites the plaintext block to declare a fake length matching
    a tampered suspect sequence. sequence_length_match passes (declared==actual)
    so the basic length check is fooled, but the ledger still has the original
    length, so ledger_length_match=False catches it."""
    fwd, rev = primers
    stamped = decoder.stamp(
        sequence=GFP_LIKE,
        synthesizer_id=51,
        run_counter=0,
        machine_state="idle",
        firmware_version="v1.0",
        primer_fwd=fwd,
        primer_rev=rev,
        ledger=fresh_ledger,
    )
    # Forge a plaintext block that declares a fake length, keeping the
    # original increment, synth_id, run_counter so the rest of the barcode
    # decodes self-consistently.
    fake_length = 9999
    forged_fields = encoder.Fields(
        increment=stamped.encoding.fields.increment,
        synthesizer_id=stamped.encoding.fields.synthesizer_id,
        run_counter=stamped.encoding.fields.run_counter,
        sequence_length=fake_length,
        mech_hash=stamped.encoding.fields.mech_hash,
        chain_hash=stamped.encoding.fields.chain_hash,
        h_sig=stamped.encoding.fields.h_sig,
        seq_hash=stamped.encoding.fields.seq_hash,
        landmarks=stamped.encoding.fields.landmarks,
    )
    fake_pt = encoder.encode_plaintext_block(forged_fields)
    pt_start = len(fwd) + layout.landmarks_after_fwd_primer
    pt_end = pt_start + layout.plaintext_len
    forged = stamped.barcode[:pt_start] + fake_pt + stamped.barcode[pt_end:]
    suspect = "A" * fake_length
    result = decoder.verify(forged, suspect, layout=layout, ledger=fresh_ledger)
    assert result.declared_sequence_length == fake_length
    assert result.sequence_length_match is True       # basic check is fooled
    assert result.ledger_sequence_length == len(GFP_LIKE)
    assert result.ledger_length_match is False        # but ledger reveals the lie
    assert "plaintext length tampered" in result.summarize()


def test_scenario3_transplant_uses_ledger_landmark_hits_for_site_forensics(
    primers, layout, fresh_ledger,
):
    """Same transplant attack, but verify is called WITHOUT
    `encoded_landmark_hits` — the realistic deployment path. The verifier
    pulls the encoder-side hits from the ledger entry that was posted at
    stamp time, so site_matches_stored is populated and we get site-level
    forensics rather than the 2-bit-feature-only fallback."""
    fwd, rev = primers
    stamped = decoder.stamp(
        sequence=GFP_LIKE,
        synthesizer_id=33,
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
        # NOTE: no encoded_landmark_hits — must come from the ledger
    )
    assert result.signature_valid
    assert not result.sequence_match
    # site-level comparison was actually performed: at least one slot anchored
    # to a different site in the suspect sequence than at encode time, AND we
    # can tell because site_matches_stored is False (not None / unevaluated).
    site_decisions = [lv.site_matches_stored for lv in result.landmarks]
    assert any(d is False for d in site_decisions), (
        "expected at least one slot where the stored site differs from the "
        "suspect-sequence site (transplant signature); got "
        f"{site_decisions}"
    )


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


def test_verify_signature_valid_is_independent_of_ledger(primers, layout, fresh_ledger):
    """After the H_sig binding fix, signature_valid is computed by reconstructing
    the (mech || seq || chain) payload from the barcode itself — the ledger
    lookup is a separate forensic axis."""
    fwd, rev = primers
    stamped = decoder.stamp(
        sequence=GFP_LIKE,
        synthesizer_id=11,
        run_counter=0,
        machine_state="idle",
        firmware_version="v1.0",
        primer_fwd=fwd,
        primer_rev=rev,
        ledger=fresh_ledger,
    )
    empty_ledger = ledger_mod.Ledger()
    result = decoder.verify(
        barcode_dna=stamped.barcode,
        suspect_sequence=GFP_LIKE,
        layout=layout,
        ledger=empty_ledger,
        encoded_landmark_hits=stamped.landmark_hits,
    )
    assert not result.ledger_hit
    assert result.signature_valid, "signature_valid must hold even without a ledger entry"


def test_verify_detects_tampered_seq_hash_via_signature(primers, layout, fresh_ledger):
    """If an attacker rewrites bases in the seq_hash region (which is not
    Hamming-protected, so a single-base change flips the recovered seq_hash),
    the verifier reconstructs (mech || seq || chain) with the new seq, signs,
    and gets a different H_sig — so signature_valid becomes False. The OLD
    'stored sig hashes to stored H_sig' check would have missed this."""
    fwd, rev = primers
    stamped = decoder.stamp(
        sequence=GFP_LIKE,
        synthesizer_id=21,
        run_counter=0,
        machine_state="idle",
        firmware_version="v1.0",
        primer_fwd=fwd,
        primer_rev=rev,
        ledger=fresh_ledger,
    )
    # locate the seq_hash bases inside the body and flip a few of them
    seq_hash_start = (
        len(fwd)
        + layout.landmarks_after_fwd_primer
        + layout.plaintext_len
        + layout.landmarks_between_plaintext_and_dict
        + layout.dict_protected_len
        + layout.landmarks_between_dict_and_seqhash
    )
    bases = list(stamped.barcode)
    swap = {"A": "C", "C": "A", "G": "T", "T": "G"}
    for i in range(layout.seq_hash_len):
        bases[seq_hash_start + i] = swap[bases[seq_hash_start + i]]
    tampered = "".join(bases)

    result = decoder.verify(
        barcode_dna=tampered,
        suspect_sequence=GFP_LIKE,
        layout=layout,
        ledger=fresh_ledger,
        encoded_landmark_hits=stamped.landmark_hits,
    )
    assert not result.signature_valid, (
        "tampering with seq_hash bases must invalidate H_sig under the new "
        "binding check"
    )


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
