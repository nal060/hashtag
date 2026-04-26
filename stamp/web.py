"""Streamlit web demo for STAMP.

Runs in a browser, deployable for free on Streamlit Community Cloud.
Imports the same `decoder`, `encoder`, `landmarks`, `ledger`, and `insert`
modules as the tkinter desktop GUI — only the presentation layer changes.

Launch locally:
    streamlit run streamlit_app.py
"""

from __future__ import annotations

import io
import random as _r
import tempfile
from pathlib import Path

import streamlit as st

from . import decoder, encoder, insert as insert_mod
from . import landmarks as lm_mod
from . import layout as layout_mod
from . import ledger as ledger_mod


# ---------- demo constants (duplicated from gui.py to avoid the tkinter import) ----------

GFP_DEMO = (
    "ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGAC"
    "GGCGACGTAAACGGCCACAAGTTCAGCGTGTCCGGCGAGGGCGAGGGCGATGCCACCTAC"
    "GGCAAGCTGACCCTGAAGTTCATCTGCACCACCGGCAAGCTGCCCGTGCCCTGGCCCACC"
    "CTCGTGACCACCCTGACCTACGGCGTGCAGTGCTTCAGCCGCTACCCCGACCACATGAAG"
    "CAGCACGACTTCTTCAAGTCCGCCATGCCCGAAGGCTACGTCCAGGAGCGCACCATCTTC"
    "TTCAAGGACGACGGCAACTACAAGACCCGCGCCGAGGTGAAGTTCGAGGGCGACACCCTG"
) * 5

MCHERRY_DEMO = (
    "ATGGTGAGCAAGGGCGAGGAGGATAACATGGCCATCATCAAGGAGTTCATGCGCTTCAAG"
    "GTGCACATGGAGGGCTCCGTGAACGGCCACGAGTTCGAGATCGAGGGCGAGGGCGAGGGC"
    "CGCCCCTACGAGGGCACCCAGACCGCCAAGCTGAAGGTGACCAAGGGTGGCCCCCTGCCC"
    "TTCGCCTGGGACATCCTGTCCCCTCAGTTCATGTACGGCTCCAAGGCCTACGTGAAGCAC"
    "CCCGCCGACATCCCCGACTACTTGAAGCTGTCCTTCCCCGAGGGCTTCAAGTGGGAGCGC"
) * 4


def _demo_primer(seed: int) -> str:
    rng = _r.Random(seed)
    while True:
        out = "".join(rng.choice("ACGT") for _ in range(20))
        if encoder.check_molecular_constraints(out)[0]:
            return out


# ---------- session state helpers ----------

_BARCODE_LEGEND = (
    ("primer",         "#cccccc", "primer (5'/3')"),
    ("plaintext",      "#a8e6e6", "plaintext: synth/run/inc/length"),
    ("dict_protected", "#d9b3ff", "dict-protected: mech/chain/H_sig"),
    ("seq_hash",       "#b3ffcc", "seq_hash"),
    ("landmark",       "#ffe680", "landmark"),
)
_COLOR_BY_TAG = {tag: c for tag, c, _ in _BARCODE_LEGEND}


def _init_state() -> None:
    if "ledger" in st.session_state:
        return
    st.session_state.ledger = ledger_mod.Ledger()
    st.session_state.fwd_primer = _demo_primer(1)
    st.session_state.rev_primer = _demo_primer(2)
    st.session_state.layout = layout_mod.BarcodeLayout(
        primer_fwd=st.session_state.fwd_primer,
        primer_rev=st.session_state.rev_primer,
    )
    st.session_state.encoded_hits = {}
    st.session_state.last_barcode = ""
    st.session_state.last_seq = ""
    st.session_state.last_increment = None
    st.session_state.last_landmark_hits = None


def _clean_seq(s: str) -> str:
    out = []
    for line in s.splitlines():
        line = line.strip()
        if line.startswith(">"):
            continue
        for ch in line.upper():
            if ch in "ACGT":
                out.append(ch)
    return "".join(out)


# ---------- HTML rendering ----------

def _segment_ranges() -> list[tuple[int, int, str]]:
    L = st.session_state.layout
    ranges: list[tuple[int, int, str]] = []
    cursor = 0
    def add(width: int, tag: str) -> None:
        nonlocal cursor
        if width:
            ranges.append((cursor, cursor + width, tag))
            cursor += width
    add(len(L.primer_fwd), "primer")
    add(L.landmarks_after_fwd_primer, "landmark")
    add(L.plaintext_len, "plaintext")
    add(L.landmarks_between_plaintext_and_dict, "landmark")
    add(L.dict_protected_len, "dict_protected")
    add(L.landmarks_between_dict_and_seqhash, "landmark")
    add(L.seq_hash_len, "seq_hash")
    add(L.landmarks_before_rev_primer, "landmark")
    add(len(L.primer_rev), "primer")
    return ranges


def _colored_barcode_html(barcode: str) -> str:
    ranges = _segment_ranges()
    if len(barcode) != ranges[-1][1]:
        return f'<pre style="font-family:monospace;font-size:13px;white-space:pre-wrap">{barcode}</pre>'
    parts = [
        f'<span style="background:{_COLOR_BY_TAG[tag]}">{barcode[s:e]}</span>'
        for s, e, tag in ranges
    ]
    return (
        '<div style="font-family:Courier New, monospace;font-size:13px;'
        'word-break:break-all;line-height:1.6">'
        + "".join(parts) + "</div>"
    )


def _legend_html() -> str:
    items = "".join(
        f'<span style="background:{c};padding:2px 6px;margin-right:8px;'
        f'border:1px solid #888;font-size:12px">{label}</span>'
        for _, c, label in _BARCODE_LEGEND
    )
    return f'<div style="margin:6px 0">{items}</div>'


def _landmark_info(barcode: str, hits) -> str:
    if not hits:
        return ""
    positions = [p for s, e, t in _segment_ranges() if t == "landmark" for p in range(s, e)]
    parts = [
        f"{slot}@{pos}→{hits[slot].upstream or '·'}"
        for slot, pos in enumerate(positions) if slot < len(hits)
    ]
    return "landmarks: " + "  ".join(parts)


def _render_verdict(result: decoder.VerifyResult) -> None:
    summary = result.summarize()
    if "CLEAN" in summary:
        st.success(summary)
    elif "NOT IN LEDGER" in summary or "INVALID" in summary or "TAMPERED" in summary or "TRANSPLANT" in summary:
        st.error(summary)
    else:
        st.warning(summary)

    text = "\n".join([
        f"  synthesizer_id      = {result.synthesizer_id}",
        f"  run_counter         = {result.run_counter}",
        f"  increment           = {result.increment}",
        f"  declared length     = {result.declared_sequence_length}",
        f"  actual length       = {result.actual_sequence_length}",
        "",
        "Hash checks:",
        f"  ledger_hit          = {result.ledger_hit}",
        f"  signature_valid     = {result.signature_valid}",
        f"  sequence_match      = {result.sequence_match}",
        f"  chain_valid         = {result.chain_valid}",
        f"  length_match        = {result.sequence_length_match}",
        f"  recovered_h_sig     = {result.recovered_h_sig.hex()}",
        f"  recovered_seq_hash  = {result.recovered_seq_hash.hex()}",
        f"  expected_seq_hash   = {result.expected_seq_hash.hex()}",
        f"  recovered_chain     = {result.recovered_chain_hash.hex()}",
        f"  expected_chain      = {result.expected_chain_hash.hex()}",
        f"  recovered_mech_hash = {result.recovered_mech_hash.hex()}",
    ])
    st.code(text, language="text")

    st.markdown("**Landmark forensic pattern**")
    rows = []
    for lv in result.landmarks:
        mark = "✓" if lv.feature_matches else "✗"
        rows.append({
            "": mark,
            "slot": lv.slot,
            "orig site": lv.stored_site or "—",
            "orig pos": lv.stored_position if lv.stored_position is not None else "—",
            "orig base": lv.stored_upstream or "—",
            "new site": lv.found_site or "—",
            "new pos": lv.found_position if lv.found_position is not None else "—",
            "new base": lv.found_upstream or "—",
            "site match": lv.site_matches_stored,
        })
    st.dataframe(rows, hide_index=True, use_container_width=True)


# ---------- actions ----------

def _do_stamp(seq: str, synth_id: int, run: int, machine_state: str, firmware: str) -> None:
    if not seq:
        st.warning("Paste or load a sequence first.")
        return
    try:
        stamped = decoder.stamp(
            sequence=seq,
            synthesizer_id=synth_id,
            run_counter=run,
            machine_state=machine_state,
            firmware_version=firmware,
            primer_fwd=st.session_state.fwd_primer,
            primer_rev=st.session_state.rev_primer,
            ledger=st.session_state.ledger,
        )
    except Exception as exc:
        st.error(f"Stamp failed: {exc}")
        return
    st.session_state.encoded_hits[(synth_id, run)] = stamped.landmark_hits
    st.session_state.last_barcode = stamped.barcode
    st.session_state.last_seq = seq
    st.session_state.last_increment = stamped.encoding.increment
    st.session_state.last_landmark_hits = stamped.landmark_hits


def _do_verify(barcode: str, seq: str) -> decoder.VerifyResult | None:
    if not barcode or not seq:
        st.warning("Both barcode and sequence required.")
        return None
    try:
        parts = st.session_state.layout.extract_full(barcode)
        pt_fields, _ = encoder.decode_plaintext_block(parts["plaintext_dna"])
        hits = st.session_state.encoded_hits.get(
            (pt_fields["synthesizer_id"], pt_fields["run_counter"])
        )
        return decoder.verify(
            barcode_dna=barcode, suspect_sequence=seq,
            layout=st.session_state.layout,
            ledger=st.session_state.ledger,
            encoded_landmark_hits=hits,
        )
    except Exception as exc:
        st.error(f"Verify failed: {exc}")
        return None


# ---------- main app ----------

def main() -> None:
    st.set_page_config(page_title="STAMP demo", layout="wide", page_icon="🧬")
    _init_state()

    st.title("🧬 STAMP — DNA Watermarking for Benchtop Synthesizers")
    st.caption(
        "Synthesis Tamper-evident Attestation and Molecular Provenance — "
        "a forensic provenance tool for benchtop DNA synthesizers."
    )

    with st.sidebar:
        st.subheader("Machine context")
        synth_id = st.number_input("synthesizer_id", value=1, min_value=0, max_value=4095)
        run_counter = st.number_input("run_counter", value=0, min_value=0, max_value=65535)
        machine_state = st.text_input("machine_state", value="idle")
        firmware = st.text_input("firmware_version", value="v1.0")
        st.divider()
        if st.button("Reset ledger", use_container_width=True):
            st.session_state.ledger = ledger_mod.Ledger()
            st.session_state.encoded_hits = {}
            st.success("Ledger cleared.")

    tab_demo, tab_plasmid, tab_scenarios = st.tabs(
        ["📝 Encode + Verify", "🧫 Plasmid insert", "🎬 Demo scenarios"]
    )

    # ---------- Tab 1: Encode + Verify ----------
    with tab_demo:
        col_e, col_v = st.columns(2)
        with col_e:
            st.subheader("Encode (synthesizer side)")
            seq_input = st.text_area(
                "Sequence (FASTA-style or raw bases)",
                value=st.session_state.get("encode_seq_input", ""),
                height=220, key="encode_seq_input",
            )
            if st.button("Stamp & post to ledger", use_container_width=True, type="primary"):
                _do_stamp(_clean_seq(seq_input), synth_id, run_counter, machine_state, firmware)

            if st.session_state.last_barcode:
                st.markdown("**Barcode (output):**")
                st.markdown(_colored_barcode_html(st.session_state.last_barcode), unsafe_allow_html=True)
                st.markdown(_legend_html(), unsafe_allow_html=True)
                info = _landmark_info(st.session_state.last_barcode, st.session_state.last_landmark_hits)
                if info:
                    st.caption(info)
                st.caption(f"increment search settled at {st.session_state.last_increment}")

        with col_v:
            st.subheader("Verify (investigator side)")
            verify_barcode = st.text_area(
                "Barcode DNA",
                value=st.session_state.last_barcode,
                height=80, key="verify_barcode_input",
            )
            verify_seq = st.text_area(
                "Suspect sequence",
                value=st.session_state.last_seq,
                height=220, key="verify_seq_input",
            )
            if st.button("Verify", use_container_width=True, type="primary"):
                result = _do_verify(_clean_seq(verify_barcode), _clean_seq(verify_seq))
                if result is not None:
                    _render_verdict(result)

    # ---------- Tab 2: Plasmid insert ----------
    with tab_plasmid:
        st.subheader("Stamp barcode into an annotated plasmid")
        upload = st.file_uploader(
            "Upload .gb / .gbk / .dna",
            type=["gb", "gbk", "genbank", "dna"],
        )
        if upload is not None:
            try:
                with tempfile.NamedTemporaryFile(
                    suffix=Path(upload.name).suffix, delete=False, mode="wb"
                ) as f:
                    f.write(upload.getvalue())
                    upload_path = f.name
                record = insert_mod.read_record(upload_path)
            except Exception as exc:
                st.error(f"Failed to read {upload.name}: {exc}")
                return
            n_features = sum(1 for f in record.features if f.type != "source")
            st.info(f"**{record.id}** — {len(record.seq)} bp · {n_features} features")

            col_p1, col_p2 = st.columns([2, 1])
            with col_p1:
                min_len = st.number_input("Candidate min run length", value=20, min_value=4, max_value=200)
                cands = insert_mod.find_candidate_sites(record, min_length=int(min_len))
                if cands:
                    options = [
                        f"pos {c.start:>6}..{c.end:<6}  {c.length:>3}× {c.base}"
                        for c in cands
                    ]
                    chosen = st.selectbox("Candidate insertion sites", options=options)
                    chosen_idx = options.index(chosen)
                    default_pos = cands[chosen_idx].start
                else:
                    st.caption("(no homopolymer runs at this min length)")
                    default_pos = 0
                position = st.number_input(
                    "Insert at position",
                    value=default_pos, min_value=0, max_value=len(record.seq),
                )
                circular = st.checkbox("Circular plasmid map", value=True)
            with col_p2:
                st.metric("Plasmid length", f"{len(record.seq)} bp")
                st.metric("Annotated features", n_features)
                st.metric("Candidate sites", len(cands))

            if st.button("Stamp + insert + render", use_container_width=True, type="primary"):
                try:
                    result = insert_mod.stamp_and_insert(
                        record=record,
                        insert_at=int(position),
                        synthesizer_id=synth_id, run_counter=run_counter,
                        machine_state=machine_state, firmware_version=firmware,
                        primer_fwd=st.session_state.fwd_primer,
                        primer_rev=st.session_state.rev_primer,
                        ledger=st.session_state.ledger,
                    )
                except Exception as exc:
                    st.error(f"Stamp failed: {exc}")
                    return

                # Write .gb to BytesIO for download
                from Bio import SeqIO
                gb_buf = io.StringIO()
                SeqIO.write(result.record, gb_buf, "genbank")
                gb_bytes = gb_buf.getvalue().encode("utf-8")

                # Render PNG to a temp file, read for inline display + download
                with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as f:
                    png_path = f.name
                insert_mod.render_to_png(result.record, png_path, circular=circular)

                st.success(
                    f"Inserted {len(result.barcode)} bp at position "
                    f"{result.insertion_position} (increment={result.increment})."
                )
                st.image(png_path, use_container_width=True)

                col_d1, col_d2 = st.columns(2)
                col_d1.download_button(
                    "Download stamped .gb",
                    data=gb_bytes,
                    file_name=f"{record.id}_stamped.gb",
                    mime="text/plain",
                    use_container_width=True,
                )
                col_d2.download_button(
                    "Download figure (PNG)",
                    data=Path(png_path).read_bytes(),
                    file_name=f"{record.id}_stamped.png",
                    mime="image/png",
                    use_container_width=True,
                )

    # ---------- Tab 3: Demo scenarios ----------
    with tab_scenarios:
        st.subheader("Spec demo scenarios")
        st.caption("Each button stamps a synthetic construct and immediately verifies it under a specific threat model.")

        c1, c2, c3 = st.columns(3)
        with c1:
            if st.button("1. Clean GFP synthesis", use_container_width=True):
                _do_stamp(GFP_DEMO, 1, 0, machine_state, firmware)
                result = _do_verify(st.session_state.last_barcode, GFP_DEMO)
                if result is not None: _render_verdict(result)
        with c2:
            if st.button("2. Modified GFP", use_container_width=True):
                _do_stamp(GFP_DEMO, 2, 0, machine_state, firmware)
                rng = _r.Random(7)
                modified = list(GFP_DEMO)
                for i in range(500):
                    modified[100 + i] = rng.choice("ACGT")
                result = _do_verify(st.session_state.last_barcode, "".join(modified))
                if result is not None: _render_verdict(result)
        with c3:
            if st.button("3. GFP→mCherry transplant", use_container_width=True):
                _do_stamp(GFP_DEMO, 3, 0, machine_state, firmware)
                result = _do_verify(st.session_state.last_barcode, MCHERRY_DEMO)
                if result is not None: _render_verdict(result)
        c4, c5, _ = st.columns(3)
        with c4:
            if st.button("4. Barcode absent (no ledger entry)", use_container_width=True):
                _do_stamp(MCHERRY_DEMO, 99, 0, machine_state, firmware)
                # wipe ledger to simulate "no attestation on file"
                st.session_state.ledger = ledger_mod.Ledger()
                st.session_state.encoded_hits = {}
                result = _do_verify(st.session_state.last_barcode, MCHERRY_DEMO)
                if result is not None: _render_verdict(result)
        with c5:
            if st.button("5. Increment-search demo", use_container_width=True):
                _do_stamp(GFP_DEMO, 0, 0, machine_state, firmware)
                if st.session_state.last_increment is not None:
                    st.info(
                        f"increment search settled at attempt {st.session_state.last_increment}; "
                        f"barcode = {len(st.session_state.last_barcode)} bases."
                    )


if __name__ == "__main__":
    main()
