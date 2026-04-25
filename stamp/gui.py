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
        self.barcode_out = scrolledtext.ScrolledText(left, height=4, wrap="char")
        self.barcode_out.pack(fill="x")

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
        self._set_text(self.barcode_out, stamped.barcode)
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
        for lv in result.landmarks:
            mark = "OK" if lv.feature_matches else "X"
            site = lv.found_site or "(none)"
            pos = lv.found_position if lv.found_position is not None else "—"
            lines.append(
                f"  [{mark}] slot {lv.slot}: stored={lv.stored_feature}  "
                f"site={site} pos={pos}  site_matches_stored={lv.site_matches_stored}"
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
        self._set_text(self.barcode_out, stamped.barcode)
        messagebox.showinfo(
            "STAMP",
            f"increment search settled at {stamped.encoding.increment} attempts. "
            f"barcode = {len(stamped.barcode)} bases.",
        )


def main() -> None:
    root = tk.Tk()
    StampApp(root)
    root.mainloop()


if __name__ == "__main__":
    main()
