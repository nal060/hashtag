# STAMP

**Synthesis Tamper-evident Attestation and Molecular Provenance**

A cryptographic DNA watermarking system for benchtop DNA synthesizers. STAMP
physically barcodes a short DNA sequence (~100–120 bases) into a non-coding
region of every synthesized construct, encoding a hardware-attested signature,
a sparse fingerprint of the synthesized sequence, and forensic landmark
features.

It is a **forensic provenance tool, not a cryptographic lock**:
- Compliant synthesis is traceable.
- Non-compliant synthesis is detectable by absence.
- Investigators get a forensic chain of custody.

The system is explicitly *not* designed to defeat determined attackers with
jailbroken synthesizers — that is a governance problem, out of scope.

## Project layout

```
stamp/
├── encoder.py     # field encoding, dictionary, nonce search, barcode assembly
├── decoder.py     # barcode extraction, decoding, Hamming correction, verification
├── landmarks.py   # restriction-site-anchored forensic landmarks (BioPython)
├── hsm_mock.py    # mock HSM signing + hash composition
├── ledger.py      # mock public ledger (in-memory / JSON)
└── gui.py         # tkinter demo GUI
```

## Quickstart

```bash
pip install -r requirements.txt
python -m stamp.gui
```

## Status

Hackathon build, 48-hour window. Not production code.
