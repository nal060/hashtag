"""STAMP encoder: field bits, DNA dictionary, sequence hash, nonce search.

Two distinct base encodings are used because of a chicken-and-egg problem:
the seed is derived from synthesizer_id + run_counter + increment, but the
seed determines the dictionary. So the seed-deriving fields (and a few other
plaintext fields useful for triage) use the **standard mapping** A=0, C=1,
G=2, T=3. Everything else uses the **seed-derived dictionary**, a random
permutation of A/C/G/T.

Hamming(31,26) is applied separately over each protected bitstream
(plaintext-protected and dict-protected) before DNA conversion.
"""

from __future__ import annotations

import hashlib
import random
from dataclasses import dataclass
from typing import Iterable

from Bio.Restriction import CommOnly
from Bio.Seq import Seq

from . import bits as bitsmod


# ---------- field schema ----------

# Plaintext-encoded (standard mapping), Hamming-protected.
# We split the 8-bit increment into two 4-bit halves at opposite ends of the
# data layout. With Hamming(31,26) blocks splitting the data at bit 26,
# placing increment_hi at the start (block 1) and increment_lo at the end
# (block 2) means *both* codeword blocks vary across the increment search.
# That's what lets the search escape primer-junction forbidden sites near
# the head and tail-homopolymers near the back of the plaintext DNA.
INCREMENT_BITS = 8
INCREMENT_HALF_BITS = 4
SYNTH_ID_BITS = 12
RUN_COUNTER_BITS = 16
SEQ_LENGTH_BITS = 16
PLAINTEXT_FIELD_ORDER = (
    ("increment_hi", INCREMENT_HALF_BITS),
    ("synthesizer_id", SYNTH_ID_BITS),
    ("run_counter", RUN_COUNTER_BITS),
    ("sequence_length", SEQ_LENGTH_BITS),
    ("increment_lo", INCREMENT_HALF_BITS),
)
PLAINTEXT_DATA_BITS = sum(w for _, w in PLAINTEXT_FIELD_ORDER)  # 52

# Dictionary-encoded, Hamming-protected.
MECH_HASH_BITS = 24
CHAIN_HASH_BITS = 24
H_SIG_BITS = 24
DICT_PROTECTED_FIELD_ORDER = (
    ("mech_hash", MECH_HASH_BITS),
    ("chain_hash", CHAIN_HASH_BITS),
    ("H_sig", H_SIG_BITS),
)
DICT_PROTECTED_DATA_BITS = sum(w for _, w in DICT_PROTECTED_FIELD_ORDER)  # 48

# Dictionary-encoded, NOT Hamming-protected.
SEQ_HASH_BITS = 24
N_LANDMARKS = 10
LANDMARK_BITS = 2  # one base per landmark
DICT_UNPROTECTED_FIELD_ORDER = (
    ("seq_hash", SEQ_HASH_BITS),
    ("landmarks", N_LANDMARKS * LANDMARK_BITS),
)
DICT_UNPROTECTED_DATA_BITS = sum(w for _, w in DICT_UNPROTECTED_FIELD_ORDER)  # 36


# ---------- molecular constraints ----------

# Source of forbidden recognition sites. Override this with any iterable of
# BioPython enzymes (e.g. AllEnzymes, or a custom RestrictionBatch). Sites
# are filtered to strict-consensus ACGT (no IUPAC ambiguity) and length
# >= _MIN_SITE_LEN to keep the blacklist tractable for the increment search.
_FORBIDDEN_ENZYMES = CommOnly
_MIN_SITE_LEN = 6

FORBIDDEN_SITES: list[str] = sorted({
    e.site.upper() for e in _FORBIDDEN_ENZYMES
    if len(e.site) >= _MIN_SITE_LEN and set(e.site.upper()) <= set("ACGT")
})


def check_molecular_constraints(dna: str) -> tuple[bool, str]:
    """Return (ok, reason). reason is empty when ok."""
    for site in FORBIDDEN_SITES:
        if site in dna:
            return False, f"forbidden site {site}"
        rc = str(Seq(site).reverse_complement())
        if rc != site and rc in dna:
            return False, f"forbidden site {rc} (rc of {site})"
    if len(dna) == 0:
        return False, "empty"
    gc = (dna.count("G") + dna.count("C")) / len(dna)
    if not (0.40 <= gc <= 0.60):
        return False, f"GC out of range ({gc:.2f})"
    for base in "ACGT":
        if base * 5 in dna:
            return False, f"homopolymer {base}*5"
    for length in range(7, 15):
        for i in range(len(dna) - length * 2):
            if dna[i:i + length] in dna[i + length:]:
                return False, f"direct repeat length {length} at pos {i}"
    return True, ""


# ---------- seed and dictionary ----------

def derive_seed(synthesizer_id: int, run_counter: int, increment: int) -> bytes:
    return hashlib.sha256(
        f"{synthesizer_id}:{run_counter}:{increment}".encode()
    ).digest()


_BASES = ("A", "C", "G", "T")


def standard_dictionary() -> dict[int, str]:
    return {0: "A", 1: "C", 2: "G", 3: "T"}


def generate_dictionary(seed: bytes) -> dict[int, str]:
    """Return a seed-derived permutation of {0,1,2,3} -> {A,C,G,T}."""
    rng = random.Random(seed)
    perm = list(_BASES)
    rng.shuffle(perm)
    return {i: b for i, b in enumerate(perm)}


def invert_dictionary(d: dict[int, str]) -> dict[str, int]:
    return {v: k for k, v in d.items()}


# ---------- scrambling ----------
#
# Without scrambling, small field values (e.g. run_counter=1) produce long
# runs of zero bits, which become long base runs (homopolymers) regardless of
# the dictionary permutation. We XOR every codeword with a fixed pseudorandom
# pattern before DNA mapping. The pattern is a public constant (no security
# value); its only job is to remove the value-content → base-content bias.

_SCRAMBLE_KEY = b"STAMP-codeword-scramble-v2"


def _scramble_pattern(length: int, key_extra: bytes = b"") -> list[int]:
    """SHA-256-derived bit stream of `length` bits, optionally re-keyed.

    `key_extra` lets callers tag a per-block-or-per-attempt suffix (e.g. the
    increment value) so the same data scrambled under different suffixes
    yields different DNA — the only way the increment search can perturb
    the plaintext block when its data bits are mostly input-fixed.
    """
    out: list[int] = []
    counter = 0
    while len(out) < length:
        chunk = hashlib.sha256(_SCRAMBLE_KEY + key_extra + counter.to_bytes(4, "big")).digest()
        for byte in chunk:
            for shift in range(7, -1, -1):
                out.append((byte >> shift) & 1)
                if len(out) >= length:
                    return out
        counter += 1
    return out


def _scramble(bit_list: list[int], key_extra: bytes = b"") -> list[int]:
    pat = _scramble_pattern(len(bit_list), key_extra)
    return [b ^ p for b, p in zip(bit_list, pat)]


def _plaintext_scramble_key(increment: int) -> bytes:
    return b"plaintext|" + increment.to_bytes(2, "big")


# ---------- bits <-> DNA ----------

def bits_to_dna(bit_list: list[int], dictionary: dict[int, str]) -> str:
    if len(bit_list) % 2 != 0:
        raise ValueError(f"need even bit count, got {len(bit_list)}")
    out: list[str] = []
    for i in range(0, len(bit_list), 2):
        idx = (bit_list[i] << 1) | bit_list[i + 1]
        out.append(dictionary[idx])
    return "".join(out)


def dna_to_bits(dna: str, dictionary: dict[int, str]) -> list[int]:
    inv = invert_dictionary(dictionary)
    out: list[int] = []
    for base in dna:
        idx = inv[base]
        out.extend([(idx >> 1) & 1, idx & 1])
    return out


# ---------- field packing ----------

def _pack_field_dict(field_order: Iterable[tuple[str, int]], values: dict[str, int]) -> list[int]:
    out: list[int] = []
    for name, width in field_order:
        out.extend(bitsmod.int_to_bits(values[name], width))
    return out


def _unpack_field_dict(field_order: Iterable[tuple[str, int]], data_bits: list[int]) -> dict[str, int]:
    out: dict[str, int] = {}
    cursor = 0
    for name, width in field_order:
        out[name] = bitsmod.bits_to_int(data_bits[cursor:cursor + width])
        cursor += width
    return out


def _pad_to_even(bit_list: list[int]) -> list[int]:
    return bit_list + [0] if len(bit_list) % 2 else bit_list


@dataclass
class Fields:
    """All per-construct field values that go into the barcode."""
    increment: int
    synthesizer_id: int
    run_counter: int
    sequence_length: int
    mech_hash: bytes        # 2 bytes
    chain_hash: bytes       # 2 bytes
    h_sig: bytes            # 2 bytes
    seq_hash: bytes         # 2 bytes
    landmarks: list[int]    # 10 ints in [0,3] — one per landmark, lower 2 bits used


def encode_plaintext_block(fields: Fields) -> str:
    """Hamming-encode the plaintext fields and convert to DNA via standard mapping.

    The scramble pad is keyed by `fields.increment` so different attempts in
    the increment search produce structurally different plaintext DNA — this
    is what lets the search escape forbidden-site collisions caused by
    input-fixed bits (synth_id / run_counter / sequence_length).
    """
    inc_hi = (fields.increment >> INCREMENT_HALF_BITS) & ((1 << INCREMENT_HALF_BITS) - 1)
    inc_lo = fields.increment & ((1 << INCREMENT_HALF_BITS) - 1)
    raw = _pack_field_dict(PLAINTEXT_FIELD_ORDER, {
        "increment_hi": inc_hi,
        "synthesizer_id": fields.synthesizer_id,
        "run_counter": fields.run_counter,
        "sequence_length": fields.sequence_length,
        "increment_lo": inc_lo,
    })
    cw = bitsmod.hamming_encode(raw)
    cw = _scramble(_pad_to_even(cw), _plaintext_scramble_key(fields.increment))
    return bits_to_dna(cw, standard_dictionary())


def decode_plaintext_block(dna: str) -> tuple[dict[str, int], list[int]]:
    """Decode the plaintext block by trying every candidate increment.

    The scramble pad depends on the increment, but the increment is itself
    inside the scrambled block. We resolve the chicken-and-egg by testing
    all 2**INCREMENT_BITS candidates and picking the one that decodes to
    itself, breaking ties by lowest total Hamming syndrome.
    """
    expected_cw_bits = bitsmod.hamming_encoded_length(PLAINTEXT_DATA_BITS)
    expected_padded = expected_cw_bits + (expected_cw_bits % 2)
    expected_dna = expected_padded // 2
    if len(dna) != expected_dna:
        raise ValueError(f"plaintext block must be {expected_dna} bases, got {len(dna)}")
    cw_scrambled = dna_to_bits(dna, standard_dictionary())

    best: tuple[int, dict[str, int], list[int]] | None = None  # (syndrome_score, fields, syndromes)
    for candidate in range(1 << INCREMENT_BITS):
        cw_padded = _scramble(cw_scrambled, _plaintext_scramble_key(candidate))
        cw = cw_padded[:expected_cw_bits]
        data, syndromes = bitsmod.hamming_decode(cw, PLAINTEXT_DATA_BITS)
        fields = _unpack_field_dict(PLAINTEXT_FIELD_ORDER, data)
        decoded_inc = (fields["increment_hi"] << INCREMENT_HALF_BITS) | fields["increment_lo"]
        if decoded_inc != candidate:
            continue
        score = sum(1 if s else 0 for s in syndromes)
        if best is None or score < best[0]:
            fields = {k: v for k, v in fields.items() if k not in ("increment_hi", "increment_lo")}
            fields["increment"] = decoded_inc
            best = (score, fields, syndromes)
            if score == 0:
                break  # clean decode — no need to keep searching
    if best is None:
        raise ValueError("no self-consistent plaintext decode (no candidate increment matched)")
    return best[1], best[2]


def encode_dict_protected_block(fields: Fields, dictionary: dict[int, str]) -> str:
    raw = (
        bitsmod.bytes_to_bits(fields.mech_hash) +
        bitsmod.bytes_to_bits(fields.chain_hash) +
        bitsmod.bytes_to_bits(fields.h_sig)
    )
    cw = bitsmod.hamming_encode(raw)
    cw = _scramble(_pad_to_even(cw))
    return bits_to_dna(cw, dictionary)


def decode_dict_protected_block(dna: str, dictionary: dict[int, str]) -> tuple[dict[str, bytes], list[int]]:
    expected_cw_bits = bitsmod.hamming_encoded_length(DICT_PROTECTED_DATA_BITS)
    expected_padded = expected_cw_bits + (expected_cw_bits % 2)
    expected_dna = expected_padded // 2
    if len(dna) != expected_dna:
        raise ValueError(f"dict-protected block must be {expected_dna} bases, got {len(dna)}")
    cw_scrambled = dna_to_bits(dna, dictionary)
    cw_padded = _scramble(cw_scrambled)
    cw = cw_padded[:expected_cw_bits]
    data, syndromes = bitsmod.hamming_decode(cw, DICT_PROTECTED_DATA_BITS)
    return {
        "mech_hash": bitsmod.bits_to_bytes(data[0:MECH_HASH_BITS]),
        "chain_hash": bitsmod.bits_to_bytes(data[MECH_HASH_BITS:MECH_HASH_BITS + CHAIN_HASH_BITS]),
        "h_sig": bitsmod.bits_to_bytes(data[MECH_HASH_BITS + CHAIN_HASH_BITS:]),
    }, syndromes


def encode_seq_hash(fields: Fields, dictionary: dict[int, str]) -> str:
    bits_list = bitsmod.bytes_to_bits(fields.seq_hash)
    bits_list = _scramble(bits_list)
    return bits_to_dna(bits_list, dictionary)


def decode_seq_hash(dna: str, dictionary: dict[int, str]) -> bytes:
    if len(dna) != SEQ_HASH_BITS // 2:
        raise ValueError(f"seq_hash must be {SEQ_HASH_BITS // 2} bases")
    bits_list = dna_to_bits(dna, dictionary)
    bits_list = _scramble(bits_list)
    return bitsmod.bits_to_bytes(bits_list)


def _landmark_shifts(seed: bytes, n: int) -> list[int]:
    """Per-position additive shifts mod 4, derived from the seed.

    We use a separate shift per landmark slot rather than the fixed XOR pad
    used elsewhere, because (a) repeated landmark values are common (the same
    upstream base shows up many times) and (b) varying these per-increment is
    the only way to make landmarks_dna content vary across the increment
    search — the landmark *values* themselves are sequence-derived and fixed.
    """
    rng = random.Random(seed + b"landmarks-shift-v1")
    return [rng.randint(0, 3) for _ in range(n)]


def encode_landmarks(landmarks: list[int], dictionary: dict[int, str], seed: bytes) -> str:
    """Each landmark is a 2-bit int. Output is one base per landmark, with
    a seed-derived per-position additive shift to break repeating-value
    homopolymers and to vary the landmark block across the increment search.
    """
    if len(landmarks) != N_LANDMARKS:
        raise ValueError(f"need {N_LANDMARKS} landmarks, got {len(landmarks)}")
    shifts = _landmark_shifts(seed, len(landmarks))
    return "".join(dictionary[(m + s) & 3] for m, s in zip(landmarks, shifts))


def decode_landmarks(dna: str, dictionary: dict[int, str], seed: bytes) -> list[int]:
    if len(dna) != N_LANDMARKS:
        raise ValueError(f"landmark block must be {N_LANDMARKS} bases")
    inv = invert_dictionary(dictionary)
    shifts = _landmark_shifts(seed, len(dna))
    return [(inv[b] - s) & 3 for b, s in zip(dna, shifts)]


# ---------- sequence hash ----------

def compute_sequence_hash(sequence: str, seed: bytes, sample_every: int = 50) -> bytes:
    """Sparse-sampled SHA-256 of the sequence, truncated to `SEQ_HASH_BITS`
    bits (24 bits per the current schema).

    Detects large rearrangements/transplants. By design misses point mutations.
    """
    n_bytes = SEQ_HASH_BITS // 8
    if len(sequence) == 0:
        return b"\x00" * n_bytes
    rng = random.Random(seed)
    n_samples = max(1, len(sequence) // sample_every)
    n_samples = min(n_samples, len(sequence))
    indexes = sorted(rng.sample(range(len(sequence)), n_samples))
    sampled = "".join(sequence[i] for i in indexes)
    return hashlib.sha256(sampled.encode()).digest()[:n_bytes]


# ---------- nonce search ----------

@dataclass
class EncodingResult:
    increment: int
    seed: bytes
    dictionary: dict[int, str]
    plaintext_dna: str
    dict_protected_dna: str
    seq_hash_dna: str
    landmarks_dna: str
    fields: Fields
    full_signature: bytes  # 32 bytes — to be posted to the ledger


def assemble_body(result: EncodingResult, layout: "BarcodeLayout | None" = None) -> str:
    """Concatenate the encoded blocks into the barcode body (no primers).

    Default layout interleaves landmarks between fields per spec. If `layout` is
    provided, it dictates the exact base-by-base composition.
    """
    if layout is not None:
        return layout.assemble(result)
    return (
        result.plaintext_dna
        + result.dict_protected_dna
        + result.seq_hash_dna
        + result.landmarks_dna
    )


def find_valid_encoding(
    sequence: str,
    synthesizer_id: int,
    run_counter: int,
    mech_hash: bytes,
    landmarks: list[int],
    primer_fwd: str,
    primer_rev: str,
    sign_fn=None,
    layout: "BarcodeLayout | None" = None,
    max_attempts: int = 256,
) -> EncodingResult:
    """Loop increment = 0..max_attempts-1; first encoding that clears
    molecular constraints wins.

    `sequence` is the user's payload DNA. `mech_hash` and `landmarks` are
    inputs because they depend on the synthesizer state and the sequence, not
    on the seed. `sign_fn` is an injectable HSM signer (default uses the mock).
    """
    if sign_fn is None:
        from . import hsm_mock as _hsm
        sign_fn = _hsm.mock_hsm_sign
    from . import hsm_mock as _hsm

    sequence_length = len(sequence)

    for increment in range(max_attempts):
        seed = derive_seed(synthesizer_id, run_counter, increment)
        dictionary = generate_dictionary(seed)

        seq_hash = compute_sequence_hash(sequence, seed)
        chain_hash = _hsm.compute_chain_hash(mech_hash, seq_hash)
        full_signature = sign_fn(mech_hash, seq_hash, chain_hash)
        h_sig = _hsm.h_sig(full_signature)

        f = Fields(
            increment=increment,
            synthesizer_id=synthesizer_id,
            run_counter=run_counter,
            sequence_length=sequence_length,
            mech_hash=mech_hash,
            chain_hash=chain_hash,
            h_sig=h_sig,
            seq_hash=seq_hash,
            landmarks=list(landmarks),
        )
        plaintext_dna = encode_plaintext_block(f)
        dict_protected_dna = encode_dict_protected_block(f, dictionary)
        seq_hash_dna = encode_seq_hash(f, dictionary)
        landmarks_dna = encode_landmarks(f.landmarks, dictionary, seed)
        result = EncodingResult(
            increment=increment,
            seed=seed,
            dictionary=dictionary,
            plaintext_dna=plaintext_dna,
            dict_protected_dna=dict_protected_dna,
            seq_hash_dna=seq_hash_dna,
            landmarks_dna=landmarks_dna,
            fields=f,
            full_signature=full_signature,
        )
        body = assemble_body(result, layout)
        full = primer_fwd + body + primer_rev
        ok, _ = check_molecular_constraints(full)
        if ok:
            return result
    raise RuntimeError(f"no valid encoding found in {max_attempts} attempts")
