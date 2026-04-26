"""Tkinter demo GUI for STAMP.

A single window with three panels:

  Left:  inputs — load sequence, set machine context + primers, encode.
  Right: verify — paste a barcode + sequence, see the multi-axis verdict.
  Bottom: scenario buttons (the 5 demo scenarios from the spec).

The GUI keeps an in-memory Ledger for the session. Each "stamp" call posts
to the ledger; "verify" reads from it.
"""

from __future__ import annotations

import random as _r
import tkinter as tk
from pathlib import Path
from tkinter import filedialog, messagebox, scrolledtext, ttk

from . import decoder, encoder, landmarks as lm_mod, layout as layout_mod, ledger as ledger_mod


# ---------- demo plasmid stubs ----------
# Real Genbank loading is straightforward but for the demo we ship two short
# synthetic constructs that resemble GFP/mCherry well enough for the spec's
# transplant scenario (different lengths, similar bases, differentiable hashes).

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


# ---------- GUI ----------

class StampApp:
    # Background colors per barcode segment. Each entry is (tag, bg, legend_label).
    _BARCODE_LEGEND: tuple[tuple[str, str, str], ...] = (
        ("primer",         "#cccccc", "primer (5'/3')"),
        ("plaintext",      "#a8e6e6", "plaintext: synth/run/inc/length"),
        ("dict_protected", "#d9b3ff", "dict-protected: mech/chain/H_sig"),
        ("seq_hash",       "#b3ffcc", "seq_hash"),
        ("landmark",       "#ffe680", "landmark"),
    )

    def __init__(self, root: tk.Tk) -> None:
        self.root = root
        root.title("STAMP demo")
        root.geometry("1100x720")

        self.ledger = ledger_mod.Ledger()
        self.fwd_primer = _demo_primer(1)
        self.rev_primer = _demo_primer(2)
        self.layout = layout_mod.BarcodeLayout(
            primer_fwd=self.fwd_primer, primer_rev=self.rev_primer,
        )
        # Keep encoded landmark hits per (synth, run) so transplant verify can
        # still do site-level forensic comparison after stamping.
        self._encoded_hits: dict[tuple[int, int], list[lm_mod.LandmarkHit]] = {}

        self._build_ui()

    # ---------- UI layout ----------

    def _build_ui(self) -> None:
        top = ttk.Frame(self.root, padding=8)
        top.pack(fill="both", expand=True)

        # Left: encode pane
        left = ttk.LabelFrame(top, text="Encode (synthesizer side)", padding=8)
        left.pack(side="left", fill="both", expand=True, padx=(0, 4))

        ttk.Label(left, text="Sequence (FASTA-style or raw bases):").pack(anchor="w")
        self.encode_seq = scrolledtext.ScrolledText(left, height=10, wrap="char")
        self.encode_seq.pack(fill="both", expand=True)

        params = ttk.Frame(left)
        params.pack(fill="x", pady=4)
        self.synth_id_var = tk.IntVar(value=1)
        self.run_counter_var = tk.IntVar(value=0)
        self.machine_state_var = tk.StringVar(value="idle")
        self.firmware_var = tk.StringVar(value="v1.0")
        for label, var, width in [
            ("synth_id", self.synth_id_var, 6),
            ("run_counter", self.run_counter_var, 6),
            ("machine_state", self.machine_state_var, 8),
            ("firmware", self.firmware_var, 6),
        ]:
            ttk.Label(params, text=label).pack(side="left", padx=(0, 2))
            ttk.Entry(params, textvariable=var, width=width).pack(side="left", padx=(0, 8))

        actions = ttk.Frame(left)
        actions.pack(fill="x", pady=2)
        ttk.Button(actions, text="Load sequence", command=self._load_sequence).pack(side="left", padx=2)
        ttk.Button(actions, text="Stamp & post to ledger", command=self._stamp).pack(side="left", padx=2)
        ttk.Button(actions, text="Copy barcode", command=self._copy_barcode).pack(side="left", padx=2)

        ttk.Label(left, text="Barcode (output):").pack(anchor="w", pady=(8, 0))
        self.barcode_out = scrolledtext.ScrolledText(
            left, height=4, wrap="char", font=("Courier New", 10),
        )
        self.barcode_out.pack(fill="x")
        self._configure_barcode_tags(self.barcode_out)

        legend = ttk.Frame(left)
        legend.pack(fill="x", pady=(4, 0))
        for tag, color, label in self._BARCODE_LEGEND:
            tk.Label(
                legend, text="  ", bg=color, relief="solid", borderwidth=1,
            ).pack(side="left", padx=(4, 2))
            tk.Label(legend, text=label, font=("Helvetica", 8)).pack(
                side="left", padx=(0, 6)
            )

        # Right: verify pane
        right = ttk.LabelFrame(top, text="Verify (investigator side)", padding=8)
        right.pack(side="left", fill="both", expand=True, padx=(4, 0))

        ttk.Label(right, text="Barcode DNA:").pack(anchor="w")
        self.verify_barcode = scrolledtext.ScrolledText(right, height=4, wrap="char")
        self.verify_barcode.pack(fill="x")

        ttk.Label(right, text="Suspect sequence:").pack(anchor="w", pady=(6, 0))
        self.verify_seq = scrolledtext.ScrolledText(right, height=10, wrap="char")
        self.verify_seq.pack(fill="both", expand=True)

        ttk.Button(right, text="Verify", command=self._verify).pack(pady=4)

        ttk.Label(right, text="Verdict:").pack(anchor="w")
        self.verify_out = scrolledtext.ScrolledText(right, height=14, wrap="word", state="disabled")
        self.verify_out.pack(fill="both", expand=True)

        # Bottom: scenarios + admin
        bottom = ttk.LabelFrame(self.root, text="Demo scenarios", padding=8)
        bottom.pack(fill="x", padx=8, pady=(0, 8))
        ttk.Button(bottom, text="1. Clean GFP synthesis", command=self._scenario1).pack(side="left", padx=2)
        ttk.Button(bottom, text="2. Modified GFP", command=self._scenario2).pack(side="left", padx=2)
        ttk.Button(bottom, text="3. GFP→mCherry transplant", command=self._scenario3).pack(side="left", padx=2)
        ttk.Button(bottom, text="4. Barcode absent", command=self._scenario4).pack(side="left", padx=2)
        ttk.Button(bottom, text="5. Increment-search demo", command=self._scenario5).pack(side="left", padx=2)
        ttk.Button(bottom, text="Reset ledger", command=self._reset_ledger).pack(side="right", padx=2)
        ttk.Button(bottom, text="Edit forbidden sites...", command=self._edit_blacklist).pack(side="right", padx=2)
        ttk.Button(bottom, text="Plasmid insert...", command=self._open_plasmid_window).pack(side="right", padx=2)

    # ---------- helpers ----------

    @staticmethod
    def _clean(seq: str) -> str:
        out_chars: list[str] = []
        for line in seq.splitlines():
            line = line.strip()
            if line.startswith(">"):
                continue
            for ch in line.upper():
                if ch in "ACGT":
                    out_chars.append(ch)
        return "".join(out_chars)

    def _configure_barcode_tags(self, widget: scrolledtext.ScrolledText) -> None:
        for tag, color, _label in self._BARCODE_LEGEND:
            widget.tag_configure(tag, background=color)

    def _segment_ranges(self) -> list[tuple[int, int, str]]:
        """[(start, end, tag), ...] for the assembled full barcode under self.layout."""
        L = self.layout
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

    def _render_barcode_highlighted(self, barcode: str) -> None:
        """Insert the barcode into self.barcode_out with per-segment background tags."""
        widget = self.barcode_out
        widget.configure(state="normal")
        widget.delete("1.0", "end")
        ranges = self._segment_ranges()
        if len(barcode) != ranges[-1][1]:
            # Length mismatch — fall back to plain text rather than mis-tagging.
            widget.insert("1.0", barcode)
            return
        for start, end, tag in ranges:
            widget.insert("end", barcode[start:end], tag)

    def _set_text(self, widget: scrolledtext.ScrolledText, text: str) -> None:
        was_disabled = str(widget.cget("state")) == "disabled"
        if was_disabled:
            widget.configure(state="normal")
        widget.delete("1.0", "end")
        widget.insert("1.0", text)
        if was_disabled:
            widget.configure(state="disabled")

    def _append_text(self, widget: scrolledtext.ScrolledText, text: str) -> None:
        widget.configure(state="normal")
        widget.insert("end", text)
        widget.configure(state="disabled")

    # ---------- actions ----------

    def _load_sequence(self) -> None:
        path = filedialog.askopenfilename(
            filetypes=[("FASTA / text", "*.fasta *.fa *.txt *.gb"), ("All files", "*.*")],
        )
        if not path:
            return
        try:
            text = Path(path).read_text()
        except Exception as exc:
            messagebox.showerror("STAMP", f"could not read {path}: {exc}")
            return
        self._set_text(self.encode_seq, self._clean(text))

    def _stamp(self) -> None:
        seq = self._clean(self.encode_seq.get("1.0", "end"))
        if not seq:
            messagebox.showwarning("STAMP", "load or paste a sequence first")
            return
        try:
            stamped = decoder.stamp(
                sequence=seq,
                synthesizer_id=int(self.synth_id_var.get()),
                run_counter=int(self.run_counter_var.get()),
                machine_state=self.machine_state_var.get(),
                firmware_version=self.firmware_var.get(),
                primer_fwd=self.fwd_primer,
                primer_rev=self.rev_primer,
                ledger=self.ledger,
            )
        except Exception as exc:
            messagebox.showerror("STAMP", f"encoding failed: {exc}")
            return
        self._encoded_hits[
            (int(self.synth_id_var.get()), int(self.run_counter_var.get()))
        ] = stamped.landmark_hits
        self._render_barcode_highlighted(stamped.barcode)
        # Pre-populate the verify pane so a one-click round-trip works
        self._set_text(self.verify_barcode, stamped.barcode)
        self._set_text(self.verify_seq, seq)
        messagebox.showinfo(
            "STAMP",
            f"Stamped: increment={stamped.encoding.increment}, "
            f"barcode {len(stamped.barcode)} bases. "
            f"Ledger entry posted for ({self.synth_id_var.get()}, {self.run_counter_var.get()}).",
        )

    def _verify(self) -> None:
        barcode = self._clean(self.verify_barcode.get("1.0", "end"))
        seq = self._clean(self.verify_seq.get("1.0", "end"))
        if not barcode or not seq:
            messagebox.showwarning("STAMP", "both barcode and sequence required")
            return
        try:
            # Best-effort: read synth/run from the barcode itself to grab any
            # encoded_hits we still have for richer landmark forensics.
            parts = self.layout.extract_full(barcode)
            pt_fields, _ = encoder.decode_plaintext_block(parts["plaintext_dna"])
            hits = self._encoded_hits.get(
                (pt_fields["synthesizer_id"], pt_fields["run_counter"])
            )
            result = decoder.verify(
                barcode_dna=barcode,
                suspect_sequence=seq,
                layout=self.layout,
                ledger=self.ledger,
                encoded_landmark_hits=hits,
            )
        except Exception as exc:
            self._set_text(self.verify_out, f"verify failed: {exc}")
            return

        self._render_verdict(result)

    def _render_verdict(self, result: decoder.VerifyResult) -> None:
        lines = [
            f"VERDICT: {result.summarize()}",
            "",
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
            "",
            "Hamming integrity:",
            f"  plaintext_clean     = {result.plaintext_hamming_clean}  (syndromes {result.plaintext_hamming_syndromes})",
            f"  dict_protected_clean= {result.dict_protected_hamming_clean}  (syndromes {result.dict_protected_hamming_syndromes})",
            "",
            "Landmark forensic pattern:",
        ]
        def _f(v, default="—"):
            return str(v) if v is not None else default
        for lv in result.landmarks:
            mark = "OK" if lv.feature_matches else "X"
            lines.append(
                f"  [{mark}] slot {lv.slot}  "
                f"orig: {_f(lv.stored_site, '(none)')} @{_f(lv.stored_position)} base={_f(lv.stored_upstream)}  |  "
                f"new: {_f(lv.found_site, '(none)')} @{_f(lv.found_position)} base={_f(lv.found_upstream)}  |  "
                f"site_match={lv.site_matches_stored}"
            )
        self._set_text(self.verify_out, "\n".join(lines))

    def _copy_barcode(self) -> None:
        barcode = self.barcode_out.get("1.0", "end").strip()
        if not barcode:
            return
        self.root.clipboard_clear()
        self.root.clipboard_append(barcode)

    def _reset_ledger(self) -> None:
        self.ledger = ledger_mod.Ledger()
        self._encoded_hits.clear()
        messagebox.showinfo("STAMP", "ledger cleared")

    def _edit_blacklist(self) -> None:
        win = tk.Toplevel(self.root)
        win.title("Forbidden restriction sites")
        win.geometry("420x420")
        ttk.Label(win, text="One site per line. Saved on close.").pack(pady=4)
        text = scrolledtext.ScrolledText(win, height=20, width=40, wrap="char")
        text.pack(fill="both", expand=True, padx=8, pady=8)
        text.insert("1.0", "\n".join(encoder.FORBIDDEN_SITES))
        def save() -> None:
            new = [l.strip().upper() for l in text.get("1.0", "end").splitlines() if l.strip()]
            encoder.FORBIDDEN_SITES.clear()
            encoder.FORBIDDEN_SITES.extend(new)
            win.destroy()
        ttk.Button(win, text="Save", command=save).pack(pady=4)

    def _open_plasmid_window(self) -> None:
        PlasmidWindow(self)

    # ---------- scenario shortcuts ----------

    def _scenario1(self) -> None:
        self.synth_id_var.set(1); self.run_counter_var.set(0)
        self._set_text(self.encode_seq, GFP_DEMO)
        self._stamp()

    def _scenario2(self) -> None:
        self.synth_id_var.set(2); self.run_counter_var.set(0)
        self._set_text(self.encode_seq, GFP_DEMO)
        self._stamp()
        # mutate a 500bp slice
        rng = _r.Random(7)
        seq = list(GFP_DEMO)
        for i in range(500):
            seq[100 + i] = rng.choice("ACGT")
        self._set_text(self.verify_seq, "".join(seq))
        self._verify()

    def _scenario3(self) -> None:
        self.synth_id_var.set(3); self.run_counter_var.set(0)
        self._set_text(self.encode_seq, GFP_DEMO)
        self._stamp()
        # paste GFP barcode onto mCherry
        self._set_text(self.verify_seq, MCHERRY_DEMO)
        self._verify()

    def _scenario4(self) -> None:
        self.synth_id_var.set(99); self.run_counter_var.set(0)
        self._set_text(self.encode_seq, MCHERRY_DEMO)
        self._stamp()
        # wipe ledger to simulate "no attestation on file"
        self.ledger = ledger_mod.Ledger()
        self._verify()

    def _scenario5(self) -> None:
        self.synth_id_var.set(0); self.run_counter_var.set(0)
        self._set_text(self.encode_seq, GFP_DEMO)
        try:
            stamped = decoder.stamp(
                sequence=GFP_DEMO,
                synthesizer_id=0, run_counter=0,
                machine_state=self.machine_state_var.get(),
                firmware_version=self.firmware_var.get(),
                primer_fwd=self.fwd_primer,
                primer_rev=self.rev_primer,
                ledger=self.ledger,
            )
        except Exception as exc:
            messagebox.showerror("STAMP", str(exc))
            return
        self._render_barcode_highlighted(stamped.barcode)
        messagebox.showinfo(
            "STAMP",
            f"increment search settled at {stamped.encoding.increment} attempts. "
            f"barcode = {len(stamped.barcode)} bases.",
        )


class PlasmidWindow:
    """Toplevel window for inserting STAMP barcodes into annotated plasmid files.

    Workflow: load a Genbank plasmid → see auto-detected candidate insertion
    sites (homopolymer runs outside any annotated feature) → pick one or type
    a position → stamp + insert → save the new .gb file → render an inline
    PNG preview (composite: full plasmid + barcode zoom; circular optional).
    """

    def __init__(self, app: "StampApp") -> None:
        self.app = app
        self.win = tk.Toplevel(app.root)
        self.win.title("STAMP — plasmid insert")
        self.win.geometry("1100x780")

        self.input_path: tk.StringVar = tk.StringVar(value="")
        self.record = None  # type: ignore[assignment]
        self.candidates: list = []
        self.position_var: tk.IntVar = tk.IntVar(value=0)
        self.candidate_min_len_var: tk.IntVar = tk.IntVar(value=20)
        self.circular_var: tk.BooleanVar = tk.BooleanVar(value=True)
        self.last_result = None  # type: ignore[assignment]
        self._preview_image = None  # keep a reference so tk doesn't GC it
        self._preview_path: str | None = None  # last rendered PNG, for resize re-render
        self._resize_after_id: str | None = None

        self._build_ui()

    def _build_ui(self) -> None:
        root = ttk.Frame(self.win, padding=8)
        root.pack(fill="both", expand=True)

        # ---- top: file load ----
        top = ttk.Frame(root)
        top.pack(fill="x")
        ttk.Button(top, text="Load .gb...", command=self._load_file).pack(side="left")
        ttk.Label(top, textvariable=self.input_path).pack(side="left", padx=8)
        self.summary_var = tk.StringVar(value="(no file loaded)")
        ttk.Label(top, textvariable=self.summary_var, foreground="#666").pack(side="left", padx=8)

        # ---- main resizable area: vertical splitter holds middle row, sequence
        # views, and preview. Drag the dividers to redistribute vertical space.
        vpw = ttk.PanedWindow(root, orient="vertical")
        vpw.pack(fill="both", expand=True, pady=(8, 0))

        # ---- middle row: horizontal splitter between candidates and parameters
        middle = ttk.PanedWindow(vpw, orient="horizontal")
        vpw.add(middle, weight=1)

        cand_frame = ttk.LabelFrame(middle, text="Candidate insertion sites", padding=6)
        middle.add(cand_frame, weight=2)
        cand_top = ttk.Frame(cand_frame)
        cand_top.pack(fill="x")
        ttk.Label(cand_top, text="min run length:").pack(side="left")
        ttk.Entry(cand_top, textvariable=self.candidate_min_len_var, width=4).pack(side="left", padx=(2, 4))
        ttk.Button(cand_top, text="Re-scan", command=self._refresh_candidates).pack(side="left")
        self.candidate_list = tk.Listbox(cand_frame, height=5, font=("Courier New", 9))
        self.candidate_list.pack(fill="both", expand=True, pady=(4, 0))
        self.candidate_list.bind("<<ListboxSelect>>", self._on_candidate_selected)

        params = ttk.LabelFrame(middle, text="Insert parameters", padding=6)
        middle.add(params, weight=1)

        def _row(label: str, var, width: int = 8) -> None:
            r = ttk.Frame(params); r.pack(fill="x", pady=1)
            ttk.Label(r, text=label, width=14).pack(side="left")
            ttk.Entry(r, textvariable=var, width=width).pack(side="left")

        _row("position", self.position_var)
        _row("synth_id", self.app.synth_id_var)
        _row("run_counter", self.app.run_counter_var)
        _row("machine_state", self.app.machine_state_var)
        _row("firmware", self.app.firmware_var)
        ttk.Checkbutton(
            params, text="circular plasmid map", variable=self.circular_var
        ).pack(anchor="w", pady=(4, 2))

        actions = ttk.Frame(params); actions.pack(fill="x", pady=(8, 0))
        ttk.Button(actions, text="Stamp & save .gb", command=self._stamp_and_save).pack(fill="x", pady=1)
        ttk.Button(actions, text="Render PNG...", command=self._render_only).pack(fill="x", pady=1)

        # ---- input/output sequence views (barcode highlighted in output) ----
        io_frame = ttk.LabelFrame(vpw, text="Sequence — input ▸ output", padding=4)
        vpw.add(io_frame, weight=1)
        self.input_seq_view = scrolledtext.ScrolledText(
            io_frame, height=3, wrap="char", font=("Courier New", 8), state="disabled",
        )
        self.input_seq_view.pack(side="left", fill="both", expand=True)
        self.output_seq_view = scrolledtext.ScrolledText(
            io_frame, height=3, wrap="char", font=("Courier New", 8), state="disabled",
        )
        self.output_seq_view.pack(side="left", fill="both", expand=True, padx=(4, 0))
        self.output_seq_view.tag_configure("barcode", background="#666666", foreground="white")

        # ---- preview ----
        preview_frame = ttk.LabelFrame(vpw, text="Preview", padding=4)
        vpw.add(preview_frame, weight=3)
        self.preview_label = tk.Label(preview_frame, background="#fafafa", anchor="center")
        self.preview_label.pack(fill="both", expand=True)
        self.preview_label.bind("<Configure>", self._on_preview_resize)

        self.status_var = tk.StringVar(value="")
        ttk.Label(root, textvariable=self.status_var, foreground="#444").pack(anchor="w")

    # ---------- file IO ----------

    def _load_file(self) -> None:
        path = filedialog.askopenfilename(
            parent=self.win,
            filetypes=[
                ("Annotated DNA", "*.gb *.gbk *.genbank *.dna *.embl"),
                ("Genbank", "*.gb *.gbk *.genbank"),
                ("SnapGene", "*.dna"),
                ("All files", "*.*"),
            ],
        )
        if not path:
            return
        try:
            from . import insert as ins_mod
            self.record = ins_mod.read_record(path)
        except Exception as exc:
            messagebox.showerror("STAMP", f"failed to read {path}: {exc}")
            return
        self.input_path.set(path)
        n_features = sum(1 for f in self.record.features if f.type != "source")
        self.summary_var.set(
            f"{self.record.id}  {len(self.record.seq)} bp  {n_features} features"
        )
        self._set_seq_view(self.input_seq_view, str(self.record.seq))
        self._set_seq_view(self.output_seq_view, "")  # cleared until next stamp
        self._refresh_candidates()

    @staticmethod
    def _set_seq_view(widget: scrolledtext.ScrolledText, text: str,
                      barcode_range: tuple[int, int] | None = None) -> None:
        widget.configure(state="normal")
        widget.delete("1.0", "end")
        widget.insert("1.0", text)
        if barcode_range is not None:
            s, e = barcode_range
            widget.tag_add("barcode", f"1.0+{s}c", f"1.0+{e}c")
        widget.configure(state="disabled")

    def _refresh_candidates(self) -> None:
        if self.record is None:
            return
        from . import insert as ins_mod
        self.candidate_list.delete(0, "end")
        self.candidates = ins_mod.find_candidate_sites(
            self.record, min_length=int(self.candidate_min_len_var.get()),
        )
        if not self.candidates:
            self.candidate_list.insert(
                "end", "(no homopolymer runs at this min length)"
            )
            return
        for c in self.candidates:
            self.candidate_list.insert(
                "end",
                f" pos {c.start:>6} .. {c.end:<6}  {c.length:>3}× {c.base}",
            )

    def _on_candidate_selected(self, event=None) -> None:
        sel = self.candidate_list.curselection()
        if not sel or not self.candidates:
            return
        idx = sel[0]
        if idx < len(self.candidates):
            self.position_var.set(self.candidates[idx].start)

    # ---------- actions ----------

    def _stamp_and_save(self) -> None:
        if self.record is None:
            messagebox.showwarning("STAMP", "load a .gb file first")
            return
        out_path = filedialog.asksaveasfilename(
            parent=self.win,
            defaultextension=".gb",
            filetypes=[("Genbank", "*.gb"), ("All files", "*.*")],
            initialfile=f"{self.record.id}_stamped.gb",
        )
        if not out_path:
            return
        from . import insert as ins_mod
        from Bio import SeqIO
        try:
            self.status_var.set("stamping...")
            self.win.update_idletasks()
            result = ins_mod.stamp_and_insert(
                record=self.record,
                insert_at=int(self.position_var.get()),
                synthesizer_id=int(self.app.synth_id_var.get()),
                run_counter=int(self.app.run_counter_var.get()),
                machine_state=self.app.machine_state_var.get(),
                firmware_version=self.app.firmware_var.get(),
                primer_fwd=self.app.fwd_primer,
                primer_rev=self.app.rev_primer,
                ledger=self.app.ledger,
            )
            SeqIO.write(result.record, out_path, "genbank")
            self.last_result = result
            # also pre-fill the main window's verify pane
            layout = self.app.layout
            full = str(result.record.seq).upper()
            fwd_pos = full.find(self.app.fwd_primer)
            if fwd_pos != -1:
                barcode = full[fwd_pos:fwd_pos + layout.full_len]
                self.app._set_text(self.app.verify_barcode, barcode)
                self.app._set_text(
                    self.app.verify_seq,
                    full[:fwd_pos] + full[fwd_pos + layout.full_len:],
                )
            self.status_var.set(
                f"wrote {out_path}  (increment={result.increment}, "
                f"inserted at {result.insertion_position})"
            )
            self._set_seq_view(
                self.output_seq_view, str(result.record.seq),
                barcode_range=(
                    result.insertion_position,
                    result.insertion_position + len(result.barcode),
                ),
            )
            # auto-render preview
            self._render_preview(result.record)
        except Exception as exc:
            self.status_var.set("")
            messagebox.showerror("STAMP", f"stamp failed: {exc}")

    def _render_only(self) -> None:
        """Render either the last stamp result, or the loaded record as-is."""
        target = None
        if self.last_result is not None:
            target = self.last_result.record
        elif self.record is not None:
            target = self.record
        if target is None:
            messagebox.showwarning("STAMP", "load a .gb file first")
            return
        out_path = filedialog.asksaveasfilename(
            parent=self.win,
            defaultextension=".png",
            filetypes=[("PNG", "*.png"), ("All files", "*.*")],
            initialfile=f"{target.id}_stamped.png",
        )
        if not out_path:
            return
        try:
            from . import insert as ins_mod
            ins_mod.render_to_png(
                target, out_path, circular=self.circular_var.get(),
            )
            self._show_preview_from_path(out_path)
            self.status_var.set(f"wrote {out_path}")
        except Exception as exc:
            messagebox.showerror("STAMP", f"render failed: {exc}")

    def _render_preview(self, record) -> None:
        """Render to a temp PNG and display inline."""
        try:
            from . import insert as ins_mod
            import tempfile
            tmp = Path(tempfile.gettempdir()) / "stamp_preview.png"
            ins_mod.render_to_png(record, tmp, circular=self.circular_var.get())
            self._show_preview_from_path(tmp)
        except Exception as exc:
            self.status_var.set(f"preview render failed: {exc}")

    def _show_preview_from_path(self, path) -> None:
        try:
            from PIL import Image, ImageTk
        except ImportError:
            self.preview_label.configure(
                text="(install Pillow to see inline previews — file written to disk OK)",
                image="",
            )
            return
        self._preview_path = path
        img = Image.open(path)
        # rescale to fill the current label size, preserving aspect ratio.
        # `resize` (unlike `thumbnail`) lets us up-scale too, so dragging the
        # preview pane bigger actually enlarges the image.
        w = self.preview_label.winfo_width() or 1000
        h = self.preview_label.winfo_height() or 380
        iw, ih = img.size
        scale = min(w / iw, h / ih)
        if scale > 0:
            img = img.resize(
                (max(int(iw * scale), 1), max(int(ih * scale), 1)),
                Image.LANCZOS,
            )
        self._preview_image = ImageTk.PhotoImage(img)
        self.preview_label.configure(image=self._preview_image, text="")

    def _on_preview_resize(self, event=None) -> None:
        """Debounced re-render so the image scales with the preview pane."""
        if not getattr(self, "_preview_path", None):
            return
        if self._resize_after_id is not None:
            self.win.after_cancel(self._resize_after_id)
        self._resize_after_id = self.win.after(
            100, lambda: self._show_preview_from_path(self._preview_path),
        )


def main() -> None:
    root = tk.Tk()
    StampApp(root)
    root.mainloop()


if __name__ == "__main__":
    main()
