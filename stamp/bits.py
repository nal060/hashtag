"""Bit-array helpers and Hamming(31,26) error correction.

Bits are represented as plain `list[int]` of 0/1 — small barcodes only, so
performance isn't a concern. All field widths are MSB-first.

Hamming(31,26) corrects any single-bit error in a 31-bit codeword. We apply
it field-by-field on the protected fields after they're packed into bits.
"""

from __future__ import annotations

from typing import Iterable


# ---------- bit packing ----------

def int_to_bits(value: int, width: int) -> list[int]:
    """MSB-first bit decomposition of `value` into exactly `width` bits."""
    if value < 0 or value >= (1 << width):
        raise ValueError(f"value {value} doesn't fit in {width} bits")
    return [(value >> (width - 1 - i)) & 1 for i in range(width)]


def bits_to_int(bits: Iterable[int]) -> int:
    """MSB-first reassembly. Inverse of int_to_bits."""
    out = 0
    for b in bits:
        out = (out << 1) | (b & 1)
    return out


def bytes_to_bits(data: bytes) -> list[int]:
    """MSB-first bit expansion of a bytes object."""
    return [(b >> (7 - i)) & 1 for b in data for i in range(8)]


def bits_to_bytes(bits: list[int]) -> bytes:
    """MSB-first bit packing into bytes. Length must be a multiple of 8."""
    if len(bits) % 8 != 0:
        raise ValueError(f"bit count {len(bits)} is not a multiple of 8")
    out = bytearray()
    for i in range(0, len(bits), 8):
        byte = 0
        for j in range(8):
            byte = (byte << 1) | (bits[i + j] & 1)
        out.append(byte)
    return bytes(out)


# ---------- Hamming(31,26) ----------
#
# Codeword layout (1-indexed positions, MSB-first when serialized):
#   p1, p2, d1, p4, d2, d3, d4, p8, d5..d11, p16, d12..d26
# Parity bit at position 2^k covers all positions whose binary index has bit k set.
# We use the 31-bit version (1..31). 26 data bits + 5 parity bits.

_PARITY_POSITIONS = (1, 2, 4, 8, 16)
_DATA_POSITIONS = tuple(p for p in range(1, 32) if p not in _PARITY_POSITIONS)


def hamming_encode_block(data_bits: list[int]) -> list[int]:
    """Encode 26 data bits → 31-bit Hamming codeword (MSB-first list)."""
    if len(data_bits) != 26:
        raise ValueError(f"need 26 data bits, got {len(data_bits)}")
    codeword = [0] * 32  # 1-indexed; index 0 unused
    for i, pos in enumerate(_DATA_POSITIONS):
        codeword[pos] = data_bits[i]
    for pp in _PARITY_POSITIONS:
        parity = 0
        for pos in range(1, 32):
            if pos == pp:
                continue
            if pos & pp:
                parity ^= codeword[pos]
        codeword[pp] = parity
    return codeword[1:]  # drop the unused index 0


def hamming_decode_block(codeword_bits: list[int]) -> tuple[list[int], int]:
    """Decode a 31-bit Hamming codeword.

    Returns (data_bits, error_position). error_position is 0 if no error
    detected, otherwise the 1-indexed position of the corrected bit.
    """
    if len(codeword_bits) != 31:
        raise ValueError(f"need 31 codeword bits, got {len(codeword_bits)}")
    cw = [0] + list(codeword_bits)  # 1-indexed
    syndrome = 0
    for pp in _PARITY_POSITIONS:
        parity = 0
        for pos in range(1, 32):
            if pos & pp:
                parity ^= cw[pos]
        if parity:
            syndrome |= pp
    if syndrome != 0:
        cw[syndrome] ^= 1
    data_bits = [cw[pos] for pos in _DATA_POSITIONS]
    return data_bits, syndrome


def hamming_encode(data_bits: list[int]) -> list[int]:
    """Encode an arbitrary-length bit array, padding to a multiple of 26.

    Returns the codeword bits. The original length is *not* embedded — callers
    must remember the original data-bit length to strip padding on decode.
    """
    out: list[int] = []
    for i in range(0, len(data_bits), 26):
        block = data_bits[i:i + 26]
        if len(block) < 26:
            block = block + [0] * (26 - len(block))
        out.extend(hamming_encode_block(block))
    return out


def hamming_decode(codeword_bits: list[int], original_data_len: int) -> tuple[list[int], list[int]]:
    """Decode a stream of 31-bit Hamming codewords.

    Returns (data_bits, syndromes) where syndromes is a list of error positions
    (0 if no error in that block). `original_data_len` strips padding off the end.
    """
    if len(codeword_bits) % 31 != 0:
        raise ValueError(f"codeword length {len(codeword_bits)} not a multiple of 31")
    out: list[int] = []
    syndromes: list[int] = []
    for i in range(0, len(codeword_bits), 31):
        block = codeword_bits[i:i + 31]
        data, syn = hamming_decode_block(block)
        out.extend(data)
        syndromes.append(syn)
    return out[:original_data_len], syndromes


def hamming_encoded_length(data_len: int) -> int:
    """Return how many codeword bits hamming_encode produces for `data_len` data bits."""
    n_blocks = (data_len + 25) // 26
    return n_blocks * 31
