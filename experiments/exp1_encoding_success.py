"""Experiment 1: Encoding success rate.

For N random plasmids of varying length, attempt STAMP encoding with
max_attempts=256 and record the increment value at which a constraint-
satisfying barcode is found. Reports the distribution and failure rate.

Output: encoding_success.csv (per-trial) and encoding_success_summary.txt.
"""

from __future__ import annotations

import csv
import random
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from stamp import encoder, hsm_mock, landmarks, layout

# ---------- experiment configuration ----------

N_TRIALS = 2000
LENGTH_RANGE = (500, 5000)   # plasmid lengths sampled uniformly in [500, 5000]
GC_RANGE = (0.40, 0.60)      # plasmid GC sampled uniformly in this range

# Constraint-clean primers found by random search; not Hamming-protected,
# but neither are real PCR primers. The asserts below catch the case where
# the forbidden-site blacklist evolves and these strings stop being clean.
PRIMER_FWD = "CAGCCTGCTCTTAAATTTTC"
PRIMER_REV = "TCGCAAGTCCACCTTGTCGC"
assert encoder.check_molecular_constraints(PRIMER_FWD)[0], (
    f"PRIMER_FWD violates constraints: {encoder.check_molecular_constraints(PRIMER_FWD)[1]}"
)
assert encoder.check_molecular_constraints(PRIMER_REV)[0], (
    f"PRIMER_REV violates constraints: {encoder.check_molecular_constraints(PRIMER_REV)[1]}"
)

OUT_DIR = Path(__file__).resolve().parent.parent / "results"
OUT_DIR.mkdir(exist_ok=True)


def _random_dna(length: int, gc: float, rng: random.Random) -> str:
    """Random DNA at target GC fraction."""
    gc_each = gc / 2.0
    at_each = (1.0 - gc) / 2.0
    return "".join(
        rng.choices("ACGT", weights=[at_each, gc_each, gc_each, at_each], k=length)
    )


def main() -> None:
    rng = random.Random(0xA17A)
    layout_obj = layout.BarcodeLayout(primer_fwd=PRIMER_FWD, primer_rev=PRIMER_REV)

    rows = []
    t0 = time.time()
    for trial in range(N_TRIALS):
        seq_len = rng.randint(*LENGTH_RANGE)
        gc = rng.uniform(*GC_RANGE)
        seq = _random_dna(seq_len, gc, rng)

        mech_hash = hsm_mock.compute_mech_hash(f"trial-{trial}", "1.0")
        hits = landmarks.find_all_landmarks(seq)
        features = landmarks.features_from_hits(hits)

        try:
            result = encoder.find_valid_encoding(
                sequence=seq,
                synthesizer_id=trial % 4096,
                run_counter=trial,
                mech_hash=mech_hash,
                landmarks=features,
                primer_fwd=PRIMER_FWD,
                primer_rev=PRIMER_REV,
                layout=layout_obj,
                max_attempts=256,
            )
            increment = result.increment
            failed = False
        except RuntimeError:
            increment = None
            failed = True

        rows.append({
            "trial": trial,
            "seq_len": seq_len,
            "gc": round(gc, 3),
            "increment": increment if increment is not None else "",
            "failed": int(failed),
        })

        if (trial + 1) % 250 == 0:
            elapsed = time.time() - t0
            print(f"  {trial + 1}/{N_TRIALS}  ({elapsed:.1f}s elapsed)")

    # ---------- write per-trial CSV ----------
    csv_path = OUT_DIR / "encoding_success.csv"
    with csv_path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["trial", "seq_len", "gc", "increment", "failed"])
        writer.writeheader()
        writer.writerows(rows)

    # ---------- summary ----------
    successes = [r["increment"] for r in rows if r["increment"] != ""]
    failures = sum(r["failed"] for r in rows)

    successes_sorted = sorted(successes)
    n_success = len(successes)
    mean = sum(successes) / n_success if n_success else float("nan")

    def pct(p: float) -> int:
        if not successes_sorted:
            return -1
        idx = min(len(successes_sorted) - 1, int(round(p / 100 * (len(successes_sorted) - 1))))
        return successes_sorted[idx]

    inc0 = sum(1 for x in successes if x == 0)
    inc_le2 = sum(1 for x in successes if x <= 2)
    inc_le5 = sum(1 for x in successes if x <= 5)

    summary_lines = [
        "STAMP encoding success rate — Experiment 1",
        "=" * 50,
        f"N trials: {N_TRIALS}",
        f"Plasmid length range: {LENGTH_RANGE[0]}-{LENGTH_RANGE[1]} bp (uniform)",
        f"GC content range: {GC_RANGE[0]:.2f}-{GC_RANGE[1]:.2f} (uniform)",
        f"Primer pair: {PRIMER_FWD} / {PRIMER_REV}",
        f"Forbidden-site list size: {len(encoder.FORBIDDEN_SITES)}",
        f"max_attempts: 256",
        "",
        "Results:",
        f"  successes: {n_success}/{N_TRIALS} ({100*n_success/N_TRIALS:.2f}%)",
        f"  failures (no valid encoding within 256 attempts): {failures}",
        "",
        "Increment-of-success distribution (over successful trials):",
        f"  increment == 0:  {inc0} ({100*inc0/n_success:.1f}%)",
        f"  increment <= 2:  {inc_le2} ({100*inc_le2/n_success:.1f}%)",
        f"  increment <= 5:  {inc_le5} ({100*inc_le5/n_success:.1f}%)",
        f"  mean:    {mean:.2f}",
        f"  median:  {pct(50)}",
        f"  p75:     {pct(75)}",
        f"  p95:     {pct(95)}",
        f"  p99:     {pct(99)}",
        f"  max:     {max(successes) if successes else 'n/a'}",
    ]
    summary = "\n".join(summary_lines)
    print()
    print(summary)
    (OUT_DIR / "encoding_success_summary.txt").write_text(summary + "\n")

    print(f"\nWrote {csv_path}")
    print(f"Wrote {OUT_DIR / 'encoding_success_summary.txt'}")


if __name__ == "__main__":
    main()
