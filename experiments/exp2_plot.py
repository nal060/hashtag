"""Plot detection power vs insertion length from experiment 2 results."""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.ticker as mtick

RESULTS_DIR = Path(__file__).resolve().parent.parent / "results"
OUT_PNG = RESULTS_DIR / "detection_power.png"


def main() -> None:
    rows = []
    with (RESULTS_DIR / "detection_power.csv").open() as f:
        for r in csv.DictReader(f):
            if r["encoding_failed"] == "1":
                continue
            rows.append(r)

    by_L: dict[int, list[dict]] = {}
    for r in rows:
        L = int(r["insertion_length"])
        by_L.setdefault(L, []).append(r)

    Ls = sorted(by_L.keys())
    seq_rate = [sum(int(t["seq_hash_flag"]) for t in by_L[L]) / len(by_L[L]) for L in Ls]
    len_rate = [sum(int(t["length_flag"]) for t in by_L[L]) / len(by_L[L]) for L in Ls]
    lmk_rate = [sum(int(t["landmark_only_flag"]) for t in by_L[L]) / len(by_L[L]) for L in Ls]

    # Plot — log-scale x except for L=0 which we drop on the log plot
    fig, ax = plt.subplots(figsize=(7.0, 4.2))
    Ls_pos = [L for L in Ls if L > 0]
    seq_pos = [seq_rate[i] for i, L in enumerate(Ls) if L > 0]
    len_pos = [len_rate[i] for i, L in enumerate(Ls) if L > 0]
    lmk_pos = [lmk_rate[i] for i, L in enumerate(Ls) if L > 0]

    ax.plot(Ls_pos, [r * 100 for r in seq_pos], marker="o", label="seq_hash flag",
            color="#2E7D32", linewidth=2)
    ax.plot(Ls_pos, [r * 100 for r in len_pos], marker="s", label="length flag",
            color="#1565C0", linewidth=2)
    ax.plot(Ls_pos, [r * 100 for r in lmk_pos], marker="^",
            label="≥1 landmark mismatch (barcode-only)",
            color="#E65100", linewidth=2)

    ax.set_xscale("log")
    ax.set_xlabel("Insertion length (bp)")
    ax.set_ylabel("Detection rate")
    ax.yaxis.set_major_formatter(mtick.PercentFormatter(decimals=0))
    ax.set_ylim(-3, 105)
    ax.set_xlim(0.8, 3500)
    ax.grid(True, which="both", alpha=0.25, linestyle=":")
    ax.set_axisbelow(True)
    ax.legend(loc="center right", framealpha=0.95, fontsize=10)

    # Title and caption-friendly text
    ax.set_title(
        "Detection power vs insertion length\n"
        f"N = {len(by_L[Ls_pos[0]])} trials per insertion length, host = 1500–3500 bp",
        fontsize=11,
    )

    # Annotate the headline barcode-only number: detection rate at 1 kb
    if 1000 in Ls_pos:
        i_1k = Ls_pos.index(1000)
        ax.annotate(
            f"{lmk_pos[i_1k]*100:.0f}% at L = 1 kb",
            xy=(1000, lmk_pos[i_1k] * 100),
            xytext=(120, 70),
            fontsize=10,
            color="#E65100",
            arrowprops=dict(arrowstyle="->", color="#E65100", alpha=0.7),
        )

    fig.tight_layout()
    fig.savefig(OUT_PNG, dpi=200)
    print(f"wrote {OUT_PNG}")


if __name__ == "__main__":
    main()
