"""Tests for stamp.insert: candidate finding, annotation-preserving insertion."""

from __future__ import annotations

from pathlib import Path

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from stamp import insert as ins, layout as layout_mod, ledger as ledger_mod


# ---------- helpers ----------

def _toy_plasmid() -> SeqRecord:
    """Tiny synthetic plasmid: 200bp with two annotated features and a
    20-base 'G' run between them outside any feature."""
    seq = (
        "ATGAAACCATGGGCATGCATGCATGCATGCATGCATGCAT"   # 0..40
        "AAAACCATGGGCATGCATGCATGCATGCATGCATGCAT"     # 40..78
        "GGGGGGGGGGGGGGGGGGGG"                       # 78..98 — homopolymer
        "ACGTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCT"       # 98..133
        "ACGTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCT"       # 133..168
        "AAGGCCTTAACCGGTTAACCGGTTAACCGGTTAA"          # 168..201
    )
    rec = SeqRecord(Seq(seq), id="toy_plasmid", name="toy", description="synthetic")
    rec.annotations["molecule_type"] = "DNA"
    rec.features = [
        SeqFeature(FeatureLocation(0, 78, strand=1), type="CDS",
                   qualifiers={"label": ["upstream_gene"]}),
        SeqFeature(FeatureLocation(98, 201, strand=1), type="CDS",
                   qualifiers={"label": ["downstream_gene"]}),
    ]
    return rec


# ---------- find_candidate_sites ----------

def test_find_candidate_sites_locates_g_run():
    rec = _toy_plasmid()
    cands = ins.find_candidate_sites(rec, min_length=10)
    assert len(cands) == 1
    c = cands[0]
    assert c.base == "G"
    assert c.start == 78
    assert c.end == 98
    assert c.length == 20


def test_find_candidate_sites_does_not_filter_by_features():
    """Runs inside annotated features are still listed — the user picks."""
    rec = _toy_plasmid()
    rec.features[0] = SeqFeature(FeatureLocation(0, 200, strand=1), type="CDS",
                                  qualifiers={"label": ["whole_plasmid"]})
    cands = ins.find_candidate_sites(rec, min_length=10)
    assert len(cands) == 1
    assert cands[0].start == 78


def test_find_candidate_sites_respects_min_length():
    rec = _toy_plasmid()
    assert ins.find_candidate_sites(rec, min_length=21) == []
    assert len(ins.find_candidate_sites(rec, min_length=20)) == 1


# ---------- annotate_barcode_segments ----------

def test_annotate_barcode_segments_covers_full_barcode():
    layout = layout_mod.BarcodeLayout(primer_fwd="A" * 20, primer_rev="C" * 20)
    feats = ins.annotate_barcode_segments(layout)
    wrapper = feats[0]
    assert int(wrapper.location.start) == 0
    assert int(wrapper.location.end) == layout.full_len
    # all sub-features should be inside the wrapper, sum of widths == full_len
    sub_widths = sum(int(f.location.end) - int(f.location.start) for f in feats[1:])
    assert sub_widths == layout.full_len


def test_annotate_barcode_segments_carries_apeinfo_colors():
    layout = layout_mod.BarcodeLayout(primer_fwd="A" * 20, primer_rev="C" * 20)
    feats = ins.annotate_barcode_segments(layout)
    for f in feats:
        assert "ApEinfo_fwdcolor" in f.qualifiers
        assert "ApEinfo_revcolor" in f.qualifiers
        assert f.qualifiers["ApEinfo_fwdcolor"][0].startswith("#")


# ---------- insert_barcode ----------

def test_insert_shifts_downstream_features_and_preserves_upstream():
    rec = _toy_plasmid()
    barcode = "X" * 50  # synthetic; the bases don't matter for this test
    out = ins.insert_barcode(rec, position=85, barcode=barcode)

    assert len(out.seq) == len(rec.seq) + 50
    assert str(out.seq[85:135]) == "X" * 50

    upstream = next(f for f in out.features if f.qualifiers.get("label") == ["upstream_gene"])
    downstream = next(f for f in out.features if f.qualifiers.get("label") == ["downstream_gene"])
    assert int(upstream.location.start) == 0
    assert int(upstream.location.end) == 78
    assert int(downstream.location.start) == 98 + 50
    assert int(downstream.location.end) == 201 + 50


def test_insert_refuses_to_split_a_feature_by_default():
    rec = _toy_plasmid()
    with pytest.raises(ins.FeatureSplitError):
        ins.insert_barcode(rec, position=20, barcode="X" * 50)


def test_insert_extends_a_feature_when_asked():
    rec = _toy_plasmid()
    out = ins.insert_barcode(rec, position=20, barcode="X" * 50, on_split_feature="extend")
    upstream = next(f for f in out.features if f.qualifiers.get("label") == ["upstream_gene"])
    assert int(upstream.location.start) == 0
    assert int(upstream.location.end) == 78 + 50


def test_insert_at_zero_and_at_end_are_valid():
    rec = _toy_plasmid()
    out_start = ins.insert_barcode(rec, position=0, barcode="X" * 5)
    assert str(out_start.seq[:5]) == "XXXXX"
    out_end = ins.insert_barcode(rec, position=len(rec.seq), barcode="X" * 5)
    assert str(out_end.seq[-5:]) == "XXXXX"


def test_insert_rejects_out_of_range_position():
    rec = _toy_plasmid()
    with pytest.raises(ValueError):
        ins.insert_barcode(rec, position=-1, barcode="X")
    with pytest.raises(ValueError):
        ins.insert_barcode(rec, position=len(rec.seq) + 1, barcode="X")


# ---------- end-to-end stamp_and_insert + verify ----------

def test_stamp_and_insert_round_trips_with_verify():
    """Full demo flow: insert a real STAMP barcode into the plasmid,
    extract the barcode region back out, and confirm verify accepts it
    against the original (pre-insertion) sequence."""
    from stamp import decoder, gui  # gui supplies the demo primers
    from stamp import encoder as _enc  # noqa: F401 — confirm import path

    rec = _toy_plasmid()
    original_seq = str(rec.seq).upper()
    ledger = ledger_mod.Ledger()
    primer_fwd = gui._demo_primer(1)
    primer_rev = gui._demo_primer(2)

    result = ins.stamp_and_insert(
        record=rec,
        insert_at=88,  # inside the G-run, outside both features
        synthesizer_id=1,
        run_counter=0,
        primer_fwd=primer_fwd,
        primer_rev=primer_rev,
        ledger=ledger,
    )

    # the barcode is exactly the inserted block
    assert str(result.record.seq[88:88 + len(result.barcode)]) == result.barcode

    # extract barcode from the inserted record by primer match and verify
    layout = layout_mod.BarcodeLayout(primer_fwd=primer_fwd, primer_rev=primer_rev)
    full = str(result.record.seq).upper()
    # locate primer_fwd
    fwd_pos = full.find(primer_fwd)
    assert fwd_pos != -1
    barcode_extracted = full[fwd_pos:fwd_pos + layout.full_len]
    suspect_seq = full[:fwd_pos] + full[fwd_pos + layout.full_len:]
    assert suspect_seq == original_seq

    verdict = decoder.verify(
        barcode_dna=barcode_extracted,
        suspect_sequence=suspect_seq,
        layout=layout,
        ledger=ledger,
    )
    assert verdict.signature_valid
    assert verdict.sequence_match
    assert verdict.chain_valid


# ---------- Genbank IO round-trip ----------

def test_render_to_png_produces_a_file(tmp_path: Path):
    pytest.importorskip("dna_features_viewer")
    from stamp import gui
    rec = _toy_plasmid()
    result = ins.stamp_and_insert(
        record=rec,
        insert_at=88,
        synthesizer_id=1,
        run_counter=0,
        primer_fwd=gui._demo_primer(1),
        primer_rev=gui._demo_primer(2),
    )
    out = tmp_path / "stamped.png"
    ins.render_to_png(result.record, out)
    assert out.exists()
    assert out.stat().st_size > 1000  # non-trivial PNG


def test_inserted_record_writes_and_reads_back(tmp_path: Path):
    rec = _toy_plasmid()
    layout = layout_mod.BarcodeLayout(primer_fwd="A" * 20, primer_rev="C" * 20)
    barcode_features = ins.annotate_barcode_segments(layout)
    out = ins.insert_barcode(rec, position=88, barcode="N" * layout.full_len,
                             barcode_features=barcode_features)
    path = tmp_path / "out.gb"
    SeqIO.write(out, path, "genbank")
    # re-read and confirm the wrapper feature survived
    re = SeqIO.read(path, "genbank")
    assert len(re.seq) == len(out.seq)
    labels = {f.qualifiers.get("label", [""])[0] for f in re.features}
    assert "STAMP barcode" in labels
    assert "STAMP plaintext block" in labels
