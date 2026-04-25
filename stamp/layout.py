"""Physical barcode layout: places encoded blocks + primers + landmark
spacers into a single DNA string, and reverses it on the way back.

Spec layout per section "Physical Layout":

  [5' primer]
  [landscape 0..2]   3 bases (clustered near primer)
  [synthesizer_id, run_counter, increment, sequence_length]  plaintext block
  [landscape 3]
  [mech_hash, chain_hash, H_sig]  dict-protected block
  [landscape 4]
  [seq_hash]
  [landscape 5..9]   5 bases (last cluster, includes those near 3' primer)
  [3' primer]

Landmark slots are interspersed as physical spacers so transplanting any
single hash chunk requires multiple independent precision operations rather
than one cut-and-paste.
"""

from __future__ import annotations

from dataclasses import dataclass

from . import encoder


@dataclass(frozen=True)
class BarcodeLayout:
    """Concrete placement of blocks within the barcode body.

    Field order is: 3 landmarks, plaintext block, 1 landmark, dict-protected
    block, 1 landmark, seq-hash block, 5 landmarks. Total of 10 landmarks
    interspersed. The exact sizes are derived from the encoder constants.
    """
    primer_fwd: str
    primer_rev: str

    # Cluster sizes (must sum to encoder.N_LANDMARKS)
    landmarks_after_fwd_primer: int = 3
    landmarks_between_plaintext_and_dict: int = 1
    landmarks_between_dict_and_seqhash: int = 1
    landmarks_before_rev_primer: int = 5

    def __post_init__(self) -> None:
        total = (
            self.landmarks_after_fwd_primer
            + self.landmarks_between_plaintext_and_dict
            + self.landmarks_between_dict_and_seqhash
            + self.landmarks_before_rev_primer
        )
        if total != encoder.N_LANDMARKS:
            raise ValueError(
                f"landmark cluster sizes sum to {total}, expected {encoder.N_LANDMARKS}"
            )

    # ---------- assembly ----------

    def assemble(self, result: encoder.EncodingResult) -> str:
        """Return the barcode body (without primers).

        The full barcode DNA is produced by `assemble_full(result)`.
        """
        lm = result.landmarks_dna
        a = self.landmarks_after_fwd_primer
        b = a + self.landmarks_between_plaintext_and_dict
        c = b + self.landmarks_between_dict_and_seqhash
        d = c + self.landmarks_before_rev_primer
        assert d == encoder.N_LANDMARKS

        return (
            lm[0:a]
            + result.plaintext_dna
            + lm[a:b]
            + result.dict_protected_dna
            + lm[b:c]
            + result.seq_hash_dna
            + lm[c:d]
        )

    def assemble_full(self, result: encoder.EncodingResult) -> str:
        return self.primer_fwd + self.assemble(result) + self.primer_rev

    # ---------- extraction ----------

    @property
    def plaintext_len(self) -> int:
        # 31 bases for 52 data bits → 2 Hamming(31,26) blocks (62 codeword bits).
        from . import bits as _bm
        cw = _bm.hamming_encoded_length(encoder.PLAINTEXT_DATA_BITS)
        return (cw + (cw % 2)) // 2

    @property
    def dict_protected_len(self) -> int:
        from . import bits as _bm
        cw = _bm.hamming_encoded_length(encoder.DICT_PROTECTED_DATA_BITS)
        return (cw + (cw % 2)) // 2

    @property
    def seq_hash_len(self) -> int:
        return encoder.SEQ_HASH_BITS // 2

    @property
    def body_len(self) -> int:
        return (
            encoder.N_LANDMARKS
            + self.plaintext_len
            + self.dict_protected_len
            + self.seq_hash_len
        )

    @property
    def full_len(self) -> int:
        return len(self.primer_fwd) + self.body_len + len(self.primer_rev)

    def extract(self, body: str) -> dict[str, str]:
        """Split the body back into named blocks and the landmarks_dna string."""
        if len(body) != self.body_len:
            raise ValueError(
                f"body must be {self.body_len} bases, got {len(body)}"
            )
        cursor = 0
        a = self.landmarks_after_fwd_primer
        b = self.landmarks_between_plaintext_and_dict
        c = self.landmarks_between_dict_and_seqhash
        d = self.landmarks_before_rev_primer

        lm_pieces: list[str] = []

        lm_pieces.append(body[cursor:cursor + a])
        cursor += a
        plaintext_dna = body[cursor:cursor + self.plaintext_len]
        cursor += self.plaintext_len

        lm_pieces.append(body[cursor:cursor + b])
        cursor += b
        dict_protected_dna = body[cursor:cursor + self.dict_protected_len]
        cursor += self.dict_protected_len

        lm_pieces.append(body[cursor:cursor + c])
        cursor += c
        seq_hash_dna = body[cursor:cursor + self.seq_hash_len]
        cursor += self.seq_hash_len

        lm_pieces.append(body[cursor:cursor + d])
        cursor += d

        if cursor != len(body):
            raise RuntimeError(f"layout extract cursor mismatch: {cursor} != {len(body)}")

        return {
            "plaintext_dna": plaintext_dna,
            "dict_protected_dna": dict_protected_dna,
            "seq_hash_dna": seq_hash_dna,
            "landmarks_dna": "".join(lm_pieces),
        }

    def extract_full(self, full: str) -> dict[str, str]:
        """Strip primers and call `extract`."""
        if len(full) != self.full_len:
            raise ValueError(f"full barcode must be {self.full_len} bases, got {len(full)}")
        if not full.startswith(self.primer_fwd):
            raise ValueError("forward primer mismatch")
        if not full.endswith(self.primer_rev):
            raise ValueError("reverse primer mismatch")
        body = full[len(self.primer_fwd):len(full) - len(self.primer_rev)]
        return self.extract(body)
