"""Experiment 2: Detection power as a function of insertion length.

For each insertion length L, run N trials of:
  1. Generate random host plasmid.
  2. STAMP it.
  3. Insert L bp of random DNA at a random host position.
  4. Verify the (unchanged) barcode against the modified host.

For each trial, record which detection signals fire:
  * sequence-hash mismatch       (seq_hash test)
  * length-field mismatch        (plaintext length test)
  * >= 1 landmark feature mismatch  (privacy-preserving 2-bit signal,
        i.e. what's recoverable from the barcode alone, with no ledger
        landmark publication)

The "any flag" detection rate is the union of the three. The "landmark-only"
detection rate isolates the privacy-preserving variant's signal.

Output: detection_power.csv (per-trial) and detection_power_summary.txt.
"""

from __future__ import annotations

import csv
import random
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from stamp import decoder, encoder, layout, ledger

# ---------- experiment configuration ----------

INSERTION_LENGTHS = [0, 1, 5, 10, 25, 50, 100, 250, 500, 1000, 2500]
N_PER_LENGTH = 200
HOST_LENGTH_RANGE = (1500, 3500)   # host plasmid lengths

PRIMER_FWD = "CAGCCTGCTCTTAAATTTTC"
PRIMER_REV = "TCGCAAGTCCACCTTGTCGC"
assert encoder.check_molecular_constraints(PRIMER_FWD)[0], "PRIMER_FWD violates constraints"
assert encoder.check_molecular_constraints(PRIMER_REV)[0], "PRIMER_REV violates constraints"

OUT_DIR = Path(__file__).resolve().parent.parent / "results"
OUT_DIR.mkdir(exist_ok=True)


def _random_dna(length: int, rng: random.Random) -> str:
    return "".join(rng.choices("ACGT", k=length))


def _count_landmark_mismatches_priv(verdict) -> int:
    """Count landmark feature mismatches under the PRIVACY-PRESERVING signal.

    This emulates verify-without-ledger-landmark-hits: for each landmark
    slot, compare the 2-bit feature stored in the barcode against the
    feature recomputed from the suspect sequence. None (no site found) is
    treated as a mismatch — this is the worst-case interpretation that the
    barcode-only verifier sees.
    """
    n = 0
    for lv in verdict.landmarks:
        found = lv.found_position is not None
        if not found:
            n += 1
            continue
        # found_hit.feature was looked up from the suspect sequence
        # We don't have direct access to it on LandmarkVerdict, so we
        # use the upstream base. found_upstream is the suspect's 5' base.
        if lv.found_upstream is None:
            n += 1
            continue
        recomputed_feature = "ACGT".index(lv.found_upstream)
        if recomputed_feature != lv.stored_feature:
            n += 1
    return n


def main() -> None:
    rng = random.Random(0xB17E)
    layout_obj = layout.BarcodeLayout(primer_fwd=PRIMER_FWD, primer_rev=PRIMER_REV)

    rows = []
    t0 = time.time()
    trial_id = 0
    for L in INSERTION_LENGTHS:
        for j in range(N_PER_LENGTH):
            host_len = rng.randint(*HOST_LENGTH_RANGE)
            host = _random_dna(host_len, rng)

            ledg = ledger.Ledger()
            try:
                stamp_result = decoder.stamp(
                    sequence=host,
                    synthesizer_id=trial_id % 4096,
                    run_counter=trial_id,
                    machine_state="clean",
                    firmware_version="1.0",
                    primer_fwd=PRIMER_FWD,
                    primer_rev=PRIMER_REV,
                    ledger=ledg,
                )
            except RuntimeError:
                # Encoding failed; skip — counted separately
                rows.append({
                    "trial": trial_id, "insertion_length": L, "host_len": host_len,
                    "encoding_failed": 1,
                    "seq_hash_flag": "", "length_flag": "",
                    "n_landmark_mismatches_priv": "",
                    "any_flag": "", "landmark_only_flag": "",
                    "verdict": "encoding_failed",
                })
                trial_id += 1
                continue

            # Apply insertion at a random host position
            if L > 0:
                ins_pos = rng.randint(0, len(host))
                insert_dna = _random_dna(L, rng)
                modified = host[:ins_pos] + insert_dna + host[ins_pos:]
            else:
                modified = host

            verdict = decoder.verify(
                stamp_result.barcode, modified, layout_obj, ledg,
                # Force barcode-only mode by NOT passing encoded_landmark_hits.
                # Note: ledger DOES contain landmark_hits because decoder.stamp
                # posted them, so the verifier auto-uses them. To get the strict
                # privacy-preserving 2-bit signal we count manually below.
                encoded_landmark_hits=None,
            )

            seq_hash_flag = not verdict.sequence_match
            length_flag = not verdict.sequence_length_match
            n_mm_priv = _count_landmark_mismatches_priv(verdict)
            landmark_only_flag = n_mm_priv >= 1
            any_flag = seq_hash_flag or length_flag or landmark_only_flag

            rows.append({
                "trial": trial_id,
                "insertion_length": L,
                "host_len": host_len,
                "encoding_failed": 0,
                "seq_hash_flag": int(seq_hash_flag),
                "length_flag": int(length_flag),
                "n_landmark_mismatches_priv": n_mm_priv,
                "any_flag": int(any_flag),
                "landmark_only_flag": int(landmark_only_flag),
                "verdict": verdict.summarize(),
            })
            trial_id += 1

        elapsed = time.time() - t0
        print(f"  L={L:>5} bp: {N_PER_LENGTH} trials done ({elapsed:.1f}s elapsed)")

    # ---------- write per-trial CSV ----------
    csv_path = OUT_DIR / "detection_power.csv"
    with csv_path.open("w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "trial", "insertion_length", "host_len", "encoding_failed",
                "seq_hash_flag", "length_flag", "n_landmark_mismatches_priv",
                "any_flag", "landmark_only_flag", "verdict",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)

    # ---------- summary table ----------
    by_L: dict[int, list[dict]] = {L: [] for L in INSERTION_LENGTHS}
    for r in rows:
        if r["encoding_failed"] == 1:
            continue
        by_L[r["insertion_length"]].append(r)

    lines = [
        "STAMP detection power vs insertion length — Experiment 2",
        "=" * 70,
        f"N per length: {N_PER_LENGTH}   host length: {HOST_LENGTH_RANGE[0]}-{HOST_LENGTH_RANGE[1]} bp",
        "",
        f"{'L (bp)':>8}  {'n':>5}  {'seq_hash':>10}  {'length':>8}  {'>=1 lmk':>9}  {'any':>6}  {'mean #lmk':>10}",
        "-" * 70,
    ]
    for L, trials in by_L.items():
        n = len(trials)
        if n == 0:
            continue
        seq = sum(t["seq_hash_flag"] for t in trials) / n
        lng = sum(t["length_flag"] for t in trials) / n
        lmk_any = sum(t["landmark_only_flag"] for t in trials) / n
        any_ = sum(t["any_flag"] for t in trials) / n
        mean_mm = sum(t["n_landmark_mismatches_priv"] for t in trials) / n
        lines.append(
            f"{L:>8}  {n:>5}  {seq*100:>9.1f}%  {lng*100:>7.1f}%  {lmk_any*100:>8.1f}%  {any_*100:>5.1f}%  {mean_mm:>10.2f}"
        )

    lines.extend([
        "",
        "Columns:",
        "  seq_hash : % of trials where the sparse sequence-hash check fires",
        "  length   : % where the plaintext length field disagrees with suspect length",
        "  >=1 lmk  : % where >=1 of 10 landmark 2-bit features mismatches",
        "             (the BARCODE-ONLY signal — privacy-preserving variant)",
        "  any      : % where any of the three signals fires (overall detection)",
        "  mean #lmk: mean number of landmark mismatches per trial (out of 10)",
    ])
    summary = "\n".join(lines)
    print()
    print(summary)
    (OUT_DIR / "detection_power_summary.txt").write_text(summary + "\n")

    print(f"\nWrote {csv_path}")
    print(f"Wrote {OUT_DIR / 'detection_power_summary.txt'}")


if __name__ == "__main__":
    main()
