"""Verify a recovered barcode against a (possibly modified) sequence.

Verification recovers the four field groups, recomputes the seq_hash from
the suspect sequence, looks the signature up in the ledger, recomputes
chain_hash from mech_hash, and walks the landmark slots. The output is a
multi-axis verdict — which fields agree, which disagree, and a forensic
summary of what kind of tampering pattern is consistent with the result.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

from . import encoder, hsm_mock, landmarks as lm_mod, layout as layout_mod, ledger as ledger_mod


# ---------- result shape ----------

@dataclass
class LandmarkVerdict:
    slot: int
    stored_feature: int
    site_present: bool
    site_matches_stored: bool
    feature_matches: bool
    found_site: Optional[str]
    found_position: Optional[int]


@dataclass
class VerifyResult:
    """Multi-axis verdict from verify().

    Each `*_valid` boolean is one independent check. The summary string
    interprets the *pattern* — see `summarize()`.
    """
    increment: int
    synthesizer_id: int
    run_counter: int
    declared_sequence_length: int
    actual_sequence_length: int

    # field-level integrity
    plaintext_hamming_clean: bool
    dict_protected_hamming_clean: bool
    plaintext_hamming_syndromes: list[int]
    dict_protected_hamming_syndromes: list[int]

    # cryptographic / hash checks
    ledger_hit: bool
    signature_valid: bool       # H_sig in barcode == SHA256(ledger.signature)[:2]
    sequence_match: bool        # seq_hash in barcode == compute_sequence_hash(suspect, seed)
    chain_valid: bool           # chain_hash in barcode == H(mech_hash || seq_hash)
    sequence_length_match: bool                  # declared (barcode) == actual (suspect)
    # declared (barcode) == ledger sequence_length. None when the ledger
    # entry is missing or pre-dates this field — backwards compat.
    ledger_length_match: Optional[bool] = None
    ledger_sequence_length: Optional[int] = None

    # forensic
    landmarks: list[LandmarkVerdict] = field(default_factory=list)

    # raw recovered values (helpful for debugging)
    recovered_seq_hash: bytes = b""
    expected_seq_hash: bytes = b""
    recovered_chain_hash: bytes = b""
    expected_chain_hash: bytes = b""
    recovered_h_sig: bytes = b""
    recovered_mech_hash: bytes = b""

    def all_clean(self) -> bool:
        return (
            self.signature_valid
            and self.sequence_match
            and self.chain_valid
            and self.sequence_length_match
            and self.ledger_length_match is not False  # None (skipped) or True both ok
            and all(l.feature_matches for l in self.landmarks)
        )

    def summarize(self) -> str:
        """Map the pattern of checks to the spec's interpretation table."""
        n_landmark_mismatches = sum(1 for l in self.landmarks if not l.feature_matches)
        if not self.ledger_hit:
            return "BARCODE NOT IN LEDGER — no synthesis attestation"
        if not self.signature_valid:
            return "SIGNATURE INVALID — jailbroken machine or forgery"
        if not self.chain_valid:
            return "CHAIN HASH INVALID — partial hash transplant suspected"
        if self.all_clean():
            return "CLEAN — all checks pass"
        # signature is valid but other checks fail
        notes: list[str] = []
        if not self.sequence_match:
            notes.append("sequence hash mismatch")
        if not self.sequence_length_match:
            notes.append(f"length mismatch (declared {self.declared_sequence_length}, "
                         f"actual {self.actual_sequence_length})")
        if self.ledger_length_match is False:
            notes.append(f"plaintext length tampered (barcode declares "
                         f"{self.declared_sequence_length}, ledger has "
                         f"{self.ledger_sequence_length})")
        if n_landmark_mismatches:
            notes.append(f"{n_landmark_mismatches}/{len(self.landmarks)} landmark features mismatch")
        if (
            not self.sequence_match
            and not self.sequence_length_match
            and n_landmark_mismatches >= len(self.landmarks) // 2
        ):
            return "TRANSPLANT SUSPECTED — H_sig valid but " + "; ".join(notes)
        return "MODIFIED — " + "; ".join(notes)


# ---------- verify ----------

def verify(
    barcode_dna: str,
    suspect_sequence: str,
    layout: layout_mod.BarcodeLayout,
    ledger: ledger_mod.Ledger,
    expected_mech_hash: Optional[bytes] = None,
    encoded_landmark_hits: Optional[list[lm_mod.LandmarkHit]] = None,
) -> VerifyResult:
    """Verify a recovered barcode against a sequence.

    `expected_mech_hash` is the value the legitimate synthesizer would have
    written, if it's known to the verifier (e.g., from a machine registry).
    When omitted we trust the mech_hash recovered from the barcode at face
    value and only check the chain self-consistency.

    `encoded_landmark_hits` (optional) lets the verifier compare
    *which sites* the encoder found, not just the bit-feature.
    """
    parts = layout.extract_full(barcode_dna)

    # 1. plaintext block — gives us increment / synthesizer_id / run_counter / sequence_length
    plaintext_fields, pt_syndromes = encoder.decode_plaintext_block(parts["plaintext_dna"])
    increment = plaintext_fields["increment"]
    synthesizer_id = plaintext_fields["synthesizer_id"]
    run_counter = plaintext_fields["run_counter"]
    declared_length = plaintext_fields["sequence_length"]

    # 2. derive seed and dictionary
    seed = encoder.derive_seed(synthesizer_id, run_counter, increment)
    dictionary = encoder.generate_dictionary(seed)

    # 3. dict-protected block — mech / chain / H_sig
    dict_fields, dp_syndromes = encoder.decode_dict_protected_block(
        parts["dict_protected_dna"], dictionary,
    )
    recovered_mech = dict_fields["mech_hash"]
    recovered_chain = dict_fields["chain_hash"]
    recovered_h_sig = dict_fields["h_sig"]

    # 4. seq_hash
    recovered_seq_hash = encoder.decode_seq_hash(parts["seq_hash_dna"], dictionary)

    # 5. landmarks
    stored_features = encoder.decode_landmarks(parts["landmarks_dna"], dictionary, seed)

    # 6. recompute seq_hash from suspect sequence with the same seed
    expected_seq_hash = encoder.compute_sequence_hash(suspect_sequence, seed)

    # 7. recompute chain_hash from (expected_mech_hash or recovered_mech) and recovered_seq_hash
    mech_for_chain = expected_mech_hash if expected_mech_hash is not None else recovered_mech
    expected_chain = hsm_mock.compute_chain_hash(mech_for_chain, recovered_seq_hash)

    # 8. signature check — reconstruct the payload the HSM should have signed
    # from the *barcode-recovered* hashes, sign it locally, and compare its
    # truncated hash to recovered H_sig. This binds H_sig to the (mech, seq,
    # chain) triple actually present in the barcode rather than just doing a
    # tautological "stored-signature hashes to stored-H_sig" check.
    expected_signature = hsm_mock.mock_hsm_sign(
        mech_hash=recovered_mech,
        seq_hash=recovered_seq_hash,
        chain_hash=recovered_chain,
    )
    signature_valid = hsm_mock.h_sig(expected_signature) == recovered_h_sig

    # Ledger lookup is now independent of the cryptographic check — it only
    # reports whether an attestation exists on file for this (synth, run).
    ledger_entry = ledger.lookup(synthesizer_id, run_counter)
    ledger_hit = ledger_entry is not None

    # 9. landmark forensic comparison. If the caller didn't pass
    # `encoded_landmark_hits`, try to pull them from the ledger entry —
    # this is the realistic deployment path where a verifier has only
    # public-ledger access. Older entries without `landmark_hits` fall
    # through to the 2-bit-feature-only comparison.
    if encoded_landmark_hits is None and ledger_entry and "landmark_hits" in ledger_entry:
        encoded_landmark_hits = [
            lm_mod.LandmarkHit(
                slot=h["slot"], site=h["site"], position=h["position"],
                upstream=h["upstream"], feature=h["feature"],
            )
            for h in ledger_entry["landmark_hits"]
        ]
    ledger_sequence_length = ledger_entry.get("sequence_length") if ledger_entry else None
    ledger_length_match = (
        declared_length == ledger_sequence_length
        if ledger_sequence_length is not None else None
    )
    cmps = lm_mod.compare_landmarks(
        stored_features=stored_features,
        sequence=suspect_sequence,
        encoded_hits=encoded_landmark_hits,
    )
    landmark_verdicts = [
        LandmarkVerdict(
            slot=c.slot,
            stored_feature=c.stored_feature,
            site_present=c.site_present,
            site_matches_stored=c.site_matches_stored,
            feature_matches=c.feature_matches,
            found_site=c.found_hit.site,
            found_position=c.found_hit.position,
        )
        for c in cmps
    ]

    return VerifyResult(
        increment=increment,
        synthesizer_id=synthesizer_id,
        run_counter=run_counter,
        declared_sequence_length=declared_length,
        actual_sequence_length=len(suspect_sequence),
        plaintext_hamming_clean=all(s == 0 for s in pt_syndromes),
        dict_protected_hamming_clean=all(s == 0 for s in dp_syndromes),
        plaintext_hamming_syndromes=pt_syndromes,
        dict_protected_hamming_syndromes=dp_syndromes,
        ledger_hit=ledger_hit,
        signature_valid=signature_valid,
        sequence_match=(recovered_seq_hash == expected_seq_hash),
        chain_valid=(recovered_chain == expected_chain),
        sequence_length_match=(declared_length == len(suspect_sequence)),
        ledger_length_match=ledger_length_match,
        ledger_sequence_length=ledger_sequence_length,
        landmarks=landmark_verdicts,
        recovered_seq_hash=recovered_seq_hash,
        expected_seq_hash=expected_seq_hash,
        recovered_chain_hash=recovered_chain,
        expected_chain_hash=expected_chain,
        recovered_h_sig=recovered_h_sig,
        recovered_mech_hash=recovered_mech,
    )


# ---------- one-shot encode-and-attest helper ----------

@dataclass
class StampResult:
    barcode: str
    encoding: encoder.EncodingResult
    landmark_hits: list[lm_mod.LandmarkHit]


def stamp(
    sequence: str,
    synthesizer_id: int,
    run_counter: int,
    machine_state: str,
    firmware_version: str,
    primer_fwd: str,
    primer_rev: str,
    ledger: ledger_mod.Ledger,
    timestamp: int = 0,
) -> StampResult:
    """End-to-end attestation: compute landmarks + mech_hash, find a valid
    encoding, post the full signature to `ledger`, and return the assembled
    barcode plus debug info.
    """
    layout = layout_mod.BarcodeLayout(primer_fwd=primer_fwd, primer_rev=primer_rev)
    mech_hash = hsm_mock.compute_mech_hash(machine_state, firmware_version)
    hits = lm_mod.find_all_landmarks(sequence)
    features = lm_mod.features_from_hits(hits)
    enc = encoder.find_valid_encoding(
        sequence=sequence,
        synthesizer_id=synthesizer_id,
        run_counter=run_counter,
        mech_hash=mech_hash,
        landmarks=features,
        primer_fwd=primer_fwd,
        primer_rev=primer_rev,
        layout=layout,
    )
    barcode = layout.assemble_full(enc)
    ledger.post(
        synthesizer_id=synthesizer_id,
        run_counter=run_counter,
        signature=enc.full_signature,
        timestamp=timestamp,
        seq_hash=enc.fields.seq_hash,
        landmark_hits=hits,
        sequence_length=enc.fields.sequence_length,
    )
    return StampResult(barcode=barcode, encoding=enc, landmark_hits=hits)
