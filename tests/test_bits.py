"""Tests for bit-packing helpers and Hamming(31,26)."""

from __future__ import annotations

import random

from stamp import bits


def test_int_bits_round_trip():
    for v, w in [(0, 4), (5, 4), (15, 4), (0xABCD, 16), (1 << 19, 20)]:
        assert bits.bits_to_int(bits.int_to_bits(v, w)) == v


def test_int_to_bits_msb_first():
    assert bits.int_to_bits(1, 4) == [0, 0, 0, 1]
    assert bits.int_to_bits(8, 4) == [1, 0, 0, 0]
    assert bits.int_to_bits(0xA, 4) == [1, 0, 1, 0]


def test_int_overflow_raises():
    import pytest
    with pytest.raises(ValueError):
        bits.int_to_bits(16, 4)
    with pytest.raises(ValueError):
        bits.int_to_bits(-1, 4)


def test_bytes_bits_round_trip():
    data = bytes(range(32))
    assert bits.bits_to_bytes(bits.bytes_to_bits(data)) == data


def test_bytes_to_bits_msb_first():
    assert bits.bytes_to_bits(b"\x80") == [1, 0, 0, 0, 0, 0, 0, 0]
    assert bits.bytes_to_bits(b"\x01") == [0, 0, 0, 0, 0, 0, 0, 1]


def test_hamming_encode_block_length():
    data = [0] * 26
    cw = bits.hamming_encode_block(data)
    assert len(cw) == 31
    assert all(b == 0 for b in cw)  # all-zero data → all-zero codeword


def test_hamming_no_error_round_trip():
    rng = random.Random(0)
    data = [rng.randint(0, 1) for _ in range(26)]
    cw = bits.hamming_encode_block(data)
    recovered, syn = bits.hamming_decode_block(cw)
    assert recovered == data
    assert syn == 0


def test_hamming_corrects_single_bit_error():
    """Flip every position, confirm it's corrected."""
    rng = random.Random(1)
    data = [rng.randint(0, 1) for _ in range(26)]
    cw = bits.hamming_encode_block(data)
    for flip_pos in range(31):
        corrupted = list(cw)
        corrupted[flip_pos] ^= 1
        recovered, syn = bits.hamming_decode_block(corrupted)
        assert recovered == data, f"failed to correct flip at {flip_pos}"
        assert syn == flip_pos + 1, f"wrong syndrome for flip at {flip_pos}"


def test_hamming_multi_block_round_trip():
    rng = random.Random(2)
    data = [rng.randint(0, 1) for _ in range(80)]  # 4 blocks worth (with padding)
    encoded = bits.hamming_encode(data)
    assert len(encoded) == bits.hamming_encoded_length(len(data))
    decoded, syndromes = bits.hamming_decode(encoded, len(data))
    assert decoded == data
    assert all(s == 0 for s in syndromes)


def test_hamming_multi_block_corrects_one_per_block():
    rng = random.Random(3)
    data = [rng.randint(0, 1) for _ in range(60)]
    encoded = bits.hamming_encode(data)
    # flip one bit in each 31-bit block
    n_blocks = len(encoded) // 31
    for block_idx in range(n_blocks):
        encoded[block_idx * 31 + (block_idx * 7) % 31] ^= 1
    decoded, syndromes = bits.hamming_decode(encoded, len(data))
    assert decoded == data
    assert all(s != 0 for s in syndromes)


def test_hamming_encoded_length_pads_to_block():
    assert bits.hamming_encoded_length(1) == 31
    assert bits.hamming_encoded_length(26) == 31
    assert bits.hamming_encoded_length(27) == 62
    assert bits.hamming_encoded_length(52) == 62
