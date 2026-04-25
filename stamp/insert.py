"""Insert a STAMP barcode into an annotated plasmid Genbank record.

The barcode goes into a "harmless" position — by convention a long
homopolymer run that lies entirely outside any annotated feature. All
features whose start position is at or past the insertion point get
shifted forward by the barcode length so the annotation map stays correct.

The barcode itself is annotated with one wrapper feature plus per-segment
sub-features (primer, plaintext, dict-protected, seq_hash, landmark). Each
segment carries `ApEinfo_fwdcolor`/`ApEinfo_revcolor` qualifiers matching
the GUI's legend, so SnapGene / ApE / Benchling render the barcode in the
same color scheme used elsewhere in the project.

Usable as a library or as a CLI:

    python -m stamp.insert INPUT.gb OUTPUT.gb \\
        --synth-id 1 --run 0 --insert-at 1234 \\
        [--list-candidates]
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Literal, Optional

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import CompoundLocation, FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from . import decoder, ledger as ledger_mod, layout as layout_mod


# ---------- color scheme (kept in sync with the GUI legend) ----------

_SEGMENT_STYLE: dict[str, dict[str, str]] = {
    "primer":         {"color": "#cccccc", "label": "STAMP primer",         "type": "primer_bind"},
    "plaintext":      {"color": "#a8e6e6", "label": "STAMP plaintext block","type": "misc_feature"},
    "dict_protected": {"color": "#d9b3ff", "label": "STAMP dict-protected", "type": "misc_feature"},
    "seq_hash":       {"color": "#b3ffcc", "label": "STAMP seq_hash",       "type": "misc_feature"},
    "landmark":       {"color": "#ffe680", "label": "STAMP landmark",       "type": "misc_feature"},
}

_WRAPPER_STYLE = {"color": "#666666", "label": "STAMP barcode", "type": "misc_feature"}


# ---------- candidate site finder ----------

@dataclass
class CandidateSite:
    start: int       # inclusive
    end: int         # exclusive
    base: str        # the homopolymer base
    length: int      # end - start
    overlapping_features: tuple[str, ...]  # labels of features that overlap (empty if clean)


def find_candidate_sites(
    record: SeqRecord,
    min_length: int = 10,
    bases: Iterable[str] = ("G", "C", "A", "T"),
    require_outside_features: bool = True,
) -> list[CandidateSite]:
    """List homopolymer runs in the plasmid that are at least `min_length`
    bases of one of `bases`. Sites that overlap any annotated feature are
    excluded when `require_outside_features=True`.
    """
    seq = str(record.seq).upper()
    n = len(seq)
    out: list[CandidateSite] = []
    feature_intervals = [
        (int(f.location.start), int(f.location.end), _label(f))
        for f in record.features
        if f.type != "source"
    ]
    i = 0
    while i < n:
        b = seq[i]
        if b not in bases:
            i += 1
            continue
        j = i
        while j < n and seq[j] == b:
            j += 1
        run_len = j - i
        if run_len >= min_length:
            overlaps = tuple(
                lbl for s, e, lbl in feature_intervals
                if not (e <= i or s >= j)
            )
            if not (require_outside_features and overlaps):
                out.append(CandidateSite(
                    start=i, end=j, base=b, length=run_len,
                    overlapping_features=overlaps,
                ))
        i = j
    return out


def _label(feature: SeqFeature) -> str:
    return feature.qualifiers.get("label", [feature.type])[0]


# ---------- barcode segment features ----------

def annotate_barcode_segments(layout: layout_mod.BarcodeLayout) -> list[SeqFeature]:
    """Return per-segment features (positions relative to barcode start) plus
    a wrapper feature that spans the whole barcode."""
    L = layout
    plan: list[tuple[int, str, str]] = []  # (width, segment_tag, suffix_for_label)

    def add(width: int, tag: str, suffix: str = "") -> None:
        if width > 0:
            plan.append((width, tag, suffix))

    add(len(L.primer_fwd), "primer", " 5'")
    add(L.landmarks_after_fwd_primer, "landmark", "")
    add(L.plaintext_len, "plaintext", "")
    add(L.landmarks_between_plaintext_and_dict, "landmark", "")
    add(L.dict_protected_len, "dict_protected", "")
    add(L.landmarks_between_dict_and_seqhash, "landmark", "")
    add(L.seq_hash_len, "seq_hash", "")
    add(L.landmarks_before_rev_primer, "landmark", "")
    add(len(L.primer_rev), "primer", " 3'")

    features: list[SeqFeature] = []
    cursor = 0
    landmark_index = 0
    for width, tag, suffix in plan:
        style = _SEGMENT_STYLE[tag]
        if tag == "landmark":
            for k in range(width):
                features.append(_segment_feature(
                    cursor + k, cursor + k + 1, tag,
                    label_override=f"{style['label']} {landmark_index}",
                ))
                landmark_index += 1
        else:
            features.append(_segment_feature(
                cursor, cursor + width, tag,
                label_override=f"{style['label']}{suffix}" if suffix else style["label"],
            ))
        cursor += width

    total = cursor
    wrapper = SeqFeature(
        FeatureLocation(0, total, strand=1),
        type=_WRAPPER_STYLE["type"],
        qualifiers={
            "label": [_WRAPPER_STYLE["label"]],
            "note": ["STAMP synthesis attestation barcode"],
            "ApEinfo_fwdcolor": [_WRAPPER_STYLE["color"]],
            "ApEinfo_revcolor": [_WRAPPER_STYLE["color"]],
        },
    )
    return [wrapper] + features


def _segment_feature(start: int, end: int, tag: str, label_override: str) -> SeqFeature:
    style = _SEGMENT_STYLE[tag]
    return SeqFeature(
        FeatureLocation(start, end, strand=1),
        type=style["type"],
        qualifiers={
            "label": [label_override],
            "note": [f"STAMP segment: {tag}"],
            "ApEinfo_fwdcolor": [style["color"]],
            "ApEinfo_revcolor": [style["color"]],
        },
    )


# ---------- insertion ----------

class FeatureSplitError(ValueError):
    """Raised when an insertion would split an existing feature and the
    caller asked us to refuse rather than extend or split."""


def insert_barcode(
    record: SeqRecord,
    position: int,
    barcode: str,
    barcode_features: Optional[list[SeqFeature]] = None,
    on_split_feature: Literal["error", "extend"] = "error",
) -> SeqRecord:
    """Return a new SeqRecord with `barcode` inserted at `position`.

    All features that lie entirely after the insertion point are shifted by
    `len(barcode)`. Features that strictly precede it are unchanged.
    Features that *span* the insertion point are extended by `len(barcode)`
    if `on_split_feature="extend"`, otherwise raise FeatureSplitError —
    the user should pick a position outside any feature.
    """
    if not 0 <= position <= len(record.seq):
        raise ValueError(f"position {position} out of range [0, {len(record.seq)}]")

    n = len(barcode)
    new_seq = record.seq[:position] + Seq(barcode) + record.seq[position:]

    new_features: list[SeqFeature] = []
    for f in record.features:
        loc = f.location
        start, end = int(loc.start), int(loc.end)
        if end <= position:
            new_features.append(f)
        elif start >= position:
            new_features.append(_shift_feature(f, n))
        else:
            if on_split_feature == "error":
                raise FeatureSplitError(
                    f"insertion at {position} would split feature "
                    f"{_label(f)!r} (location {start}..{end}); pick a "
                    f"position outside any annotated feature, or pass "
                    f"on_split_feature='extend'"
                )
            new_features.append(_extend_feature(f, position, n))

    if barcode_features:
        for bf in barcode_features:
            new_features.append(_shift_feature(bf, position))

    new_features.sort(key=lambda f: int(f.location.start))

    annotations = dict(record.annotations)
    annotations.setdefault("molecule_type", "DNA")  # required by Genbank writer

    new_record = SeqRecord(
        seq=new_seq,
        id=record.id,
        name=record.name,
        description=record.description,
        features=new_features,
        annotations=annotations,
    )
    new_record.dbxrefs = list(record.dbxrefs)
    return new_record


def _shift_location(loc, offset: int):
    if isinstance(loc, CompoundLocation):
        return CompoundLocation(
            [_shift_location(p, offset) for p in loc.parts],
            operator=loc.operator,
        )
    return FeatureLocation(loc.start + offset, loc.end + offset, strand=loc.strand)


def _shift_feature(f: SeqFeature, offset: int) -> SeqFeature:
    return SeqFeature(
        _shift_location(f.location, offset),
        type=f.type,
        qualifiers=dict(f.qualifiers),
    )


def _extend_feature(f: SeqFeature, position: int, n: int) -> SeqFeature:
    """Extend a feature that spans `position` by `n` bases (insertion is
    treated as part of the feature)."""
    loc = f.location
    if isinstance(loc, CompoundLocation):
        new_parts = []
        for p in loc.parts:
            ps, pe = int(p.start), int(p.end)
            if pe <= position:
                new_parts.append(p)
            elif ps >= position:
                new_parts.append(FeatureLocation(p.start + n, p.end + n, strand=p.strand))
            else:
                new_parts.append(FeatureLocation(p.start, p.end + n, strand=p.strand))
        new_loc = CompoundLocation(new_parts, operator=loc.operator)
    else:
        new_loc = FeatureLocation(loc.start, loc.end + n, strand=loc.strand)
    return SeqFeature(new_loc, type=f.type, qualifiers=dict(f.qualifiers))


# ---------- end-to-end stamp + insert ----------

@dataclass
class InsertResult:
    record: SeqRecord
    barcode: str
    insertion_position: int
    increment: int


def stamp_and_insert(
    record: SeqRecord,
    insert_at: int,
    synthesizer_id: int,
    run_counter: int,
    machine_state: str = "idle",
    firmware_version: str = "v1.0",
    primer_fwd: Optional[str] = None,
    primer_rev: Optional[str] = None,
    ledger: Optional[ledger_mod.Ledger] = None,
    on_split_feature: Literal["error", "extend"] = "error",
) -> InsertResult:
    """Generate a STAMP barcode for `record.seq` and insert it at `insert_at`.

    `record` is treated as the user's payload sequence — the barcode's
    seq_hash is computed over the *pre-insertion* sequence, which matches
    what the verifier reconstructs after PCR-amplifying the barcode out.

    `primer_fwd`/`primer_rev` default to the GUI's demo primers. Pass
    manufacturer-specific primers for real use.
    """
    if primer_fwd is None or primer_rev is None:
        from .gui import _demo_primer  # type: ignore
        if primer_fwd is None:
            primer_fwd = _demo_primer(1)
        if primer_rev is None:
            primer_rev = _demo_primer(2)
    if ledger is None:
        ledger = ledger_mod.Ledger()

    stamped = decoder.stamp(
        sequence=str(record.seq).upper(),
        synthesizer_id=synthesizer_id,
        run_counter=run_counter,
        machine_state=machine_state,
        firmware_version=firmware_version,
        primer_fwd=primer_fwd,
        primer_rev=primer_rev,
        ledger=ledger,
    )

    layout = layout_mod.BarcodeLayout(primer_fwd=primer_fwd, primer_rev=primer_rev)
    barcode_features = annotate_barcode_segments(layout)

    new_record = insert_barcode(
        record=record,
        position=insert_at,
        barcode=stamped.barcode,
        barcode_features=barcode_features,
        on_split_feature=on_split_feature,
    )
    return InsertResult(
        record=new_record,
        barcode=stamped.barcode,
        insertion_position=insert_at,
        increment=stamped.encoding.increment,
    )


# ---------- CLI ----------

def _build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(prog="stamp.insert", description=__doc__.split("\n\n")[0])
    ap.add_argument("input", help="input Genbank file (annotated plasmid)")
    ap.add_argument("output", nargs="?", help="output Genbank file (omit when using --list-candidates)")
    ap.add_argument("--synth-id", type=int, default=1, help="synthesizer_id field")
    ap.add_argument("--run", type=int, default=0, help="run_counter field")
    ap.add_argument("--insert-at", type=int, help="insertion position (0-indexed)")
    ap.add_argument("--state", default="idle", help="machine_state for mech_hash")
    ap.add_argument("--firmware", default="v1.0", help="firmware_version for mech_hash")
    ap.add_argument("--ledger", help="optional path to a JSON ledger file (created if absent)")
    ap.add_argument("--primer-fwd", help="forward primer (defaults to demo primer)")
    ap.add_argument("--primer-rev", help="reverse primer (defaults to demo primer)")
    ap.add_argument("--on-split", choices=("error", "extend"), default="error",
                    help="behavior if insertion would split a feature")
    ap.add_argument("--list-candidates", action="store_true",
                    help="print homopolymer candidate sites and exit (no insertion)")
    ap.add_argument("--candidate-min-len", type=int, default=10,
                    help="minimum homopolymer run length when listing candidates")
    return ap


def main(argv: Optional[list[str]] = None) -> int:
    args = _build_parser().parse_args(argv)
    record = SeqIO.read(args.input, "genbank")

    if args.list_candidates:
        cands = find_candidate_sites(record, min_length=args.candidate_min_len)
        if not cands:
            print(f"No homopolymer runs of length >= {args.candidate_min_len} "
                  f"found outside annotated features.")
            return 0
        print(f"{len(cands)} candidate site(s) in {record.id}:")
        for c in cands:
            print(f"  pos {c.start}..{c.end}  {c.length}× {c.base}")
        return 0

    if args.output is None or args.insert_at is None:
        print("error: --insert-at and OUTPUT are required unless --list-candidates is given")
        return 2

    ledger = ledger_mod.Ledger(Path(args.ledger)) if args.ledger else ledger_mod.Ledger()

    result = stamp_and_insert(
        record=record,
        insert_at=args.insert_at,
        synthesizer_id=args.synth_id,
        run_counter=args.run,
        machine_state=args.state,
        firmware_version=args.firmware,
        primer_fwd=args.primer_fwd,
        primer_rev=args.primer_rev,
        ledger=ledger,
        on_split_feature=args.on_split,
    )
    SeqIO.write(result.record, args.output, "genbank")
    print(f"wrote {args.output}: inserted {len(result.barcode)} bp at position "
          f"{result.insertion_position} (increment={result.increment})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
