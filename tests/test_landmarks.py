"""Tests for stamp.landmarks: priority lists, finder, comparison."""

from __future__ import annotations

from stamp import landmarks


def test_priority_lists_partition_canonical_6mers():
    piles = landmarks.generate_priority_lists()
    assert len(piles) == landmarks.N_LANDMARKS
    total = sum(len(p) for p in piles)
    # Number of canonical 6-mers: 4^6 = 4096, but pairs of strand-equivalent
    # 6-mers collapse. Should be > 2000.
    assert total > 2000
    # No 6-mer appears in more than one pile.
    seen: set[str] = set()
    for pile in piles:
        for mer in pile:
            assert mer not in seen, f"{mer} duplicated"
            seen.add(mer)


def test_priority_lists_pin_default_sites():
    piles = landmarks.generate_priority_lists()
    canonical_pinned = [landmarks._canonical(s) for s in landmarks.DEFAULT_PRIORITY_SITES]
    for i, expected_top in enumerate(canonical_pinned):
        assert piles[i][0] == expected_top, (
            f"slot {i}: expected {expected_top} at top, got {piles[i][0]}"
        )


def test_priority_lists_deterministic():
    a = landmarks.generate_priority_lists()
    b = landmarks.generate_priority_lists()
    assert a == b


# ---------- find_landmark ----------

def test_find_landmark_returns_first_priority_when_present():
    # build a sequence that contains EcoRI site (GAATTC) at position 10
    upstream = "ACGTACGTAC"  # the base immediately 5' of GAATTC is 'C'
    seq = upstream + "GAATTC" + "ACGTACGT" * 5
    piles = landmarks.generate_priority_lists()
    hit = landmarks.find_landmark(seq, piles[0], slot=0)
    assert hit.site == "GAATTC"  # already canonical
    assert hit.position == 10
    assert hit.upstream == "C"
    assert hit.feature == landmarks._BASE_TO_INT["C"]


def test_find_landmark_falls_through_when_first_absent():
    # No EcoRI site, but include the second-priority for slot 0's pile
    piles = landmarks.generate_priority_lists()
    second = piles[0][1]
    seq = "ACGT" + second + "ACGTACGT" * 5
    hit = landmarks.find_landmark(seq, piles[0], slot=0)
    assert hit.site == second
    # If second site happens to ALSO not be present, we'd fall through; but our
    # construction guarantees it.


def test_find_landmark_returns_none_when_sequence_too_short():
    """If the sequence is shorter than 6 bases, no 6-mer can be found."""
    piles = landmarks.generate_priority_lists()
    hit = landmarks.find_landmark("ACGT", piles[0], slot=0)
    assert hit.site is None
    assert hit.feature is None


def test_find_landmark_tiebreak_by_downstream():
    # Two GAATTC sites; choose the one with alphabetically smaller downstream
    seq = "ACGTACGTACGAATTCZZZZZZZZZZACGAATTCAAAAAAAAAAAA".replace("Z", "G")
    # Now positions of GAATTC: find them
    expected_positions = []
    p = 0
    while True:
        idx = seq.find("GAATTC", p)
        if idx == -1:
            break
        expected_positions.append(idx)
        p = idx + 1
    assert len(expected_positions) >= 2
    piles = landmarks.generate_priority_lists()
    hit = landmarks.find_landmark(seq, piles[0], slot=0)
    # downstream-sequence tiebreak picks the one with smaller seq[p+6:p+16]
    expected = min(
        expected_positions,
        key=lambda p: seq[p + 6:p + 16] if p + 16 <= len(seq) else seq[p + 6:],
    )
    assert hit.position == expected


def test_find_all_landmarks_returns_one_per_slot():
    seq = "GCATGCATGCATGCATGCAT" * 20
    hits = landmarks.find_all_landmarks(seq)
    assert len(hits) == landmarks.N_LANDMARKS
    for i, h in enumerate(hits):
        assert h.slot == i


def test_features_from_hits_extracts_2bit_values():
    seq = "GCATGCATGCATGCATGCAT" * 20
    hits = landmarks.find_all_landmarks(seq)
    feats = landmarks.features_from_hits(hits)
    assert len(feats) == landmarks.N_LANDMARKS
    for f in feats:
        assert 0 <= f <= 3


# ---------- compare_landmarks ----------

def test_compare_landmarks_identical_sequence_all_match():
    seq = "GCATGCATGCATGCATGCAT" * 20
    hits = landmarks.find_all_landmarks(seq)
    feats = landmarks.features_from_hits(hits)
    cmp = landmarks.compare_landmarks(feats, seq, encoded_hits=hits)
    assert all(c.feature_matches for c in cmp)
    assert all(c.site_matches_stored for c in cmp)


def test_compare_landmarks_modified_sequence_some_mismatch():
    """Replace a chunk of the sequence and at least one landmark should change."""
    seq = "GCATGCATGCATGCATGCAT" * 50
    feats = landmarks.features_from_hits(landmarks.find_all_landmarks(seq))
    modified = seq[:200] + "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT" + seq[240:]
    cmp = landmarks.compare_landmarks(feats, modified)
    n_mismatch = sum(1 for c in cmp if not c.feature_matches)
    assert n_mismatch >= 1
