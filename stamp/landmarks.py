"""Restriction-site-anchored forensic landmark features.

Landmarks are NOT cryptographically secure. Their value is forensic: they
make hash-transplantation require precision molecular biology on multiple
independent targets (each landmark is anchored to a different restriction
site in the construct), and the *pattern* of which landmarks survive a
modification tells investigators what kind of modification occurred.

Algorithm summary:
- Generate 10 priority lists of canonical 6-mers (one list per landmark slot).
- Slots 0-3 are bubble-sorted with one anchored restriction enzyme site each
  (EcoRI / BamHI / HindIII / XhoI). Slots 4-9 keep the natural round-robin
  order, so the 10 landmarks don't all collapse to a multiple-cloning-site
  cluster on plasmids that happen to carry one.
- For each landmark, scan the construct for the highest-priority site that
  is present; tiebreak alphabetically by downstream sequence.
- Read the 2 upstream bases of the chosen site; the first of those 2 bases
  (2 bits) is the landmark feature.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import product
from typing import Optional

from Bio.Seq import Seq


N_LANDMARKS = 10

# Sites anchored at the front of slots 0..len(DEFAULT_PRIORITY_SITES)-1.
# Slots beyond this length use the natural round-robin order so all 10
# landmarks don't collapse to a single MCS cluster.
DEFAULT_PRIORITY_SITES: tuple[str, ...] = (
    "GAATTC",  # EcoRI  → slot 0
    "GGATCC",  # BamHI  → slot 1
    "AAGCTT",  # HindIII → slot 2
    "CTCGAG",  # XhoI   → slot 3
)

_BASE_TO_INT = {"A": 0, "C": 1, "G": 2, "T": 3}


def _canonical(seq: str) -> str:
    """Lexicographically smaller of seq and its reverse-complement."""
    rc = str(Seq(seq).reverse_complement())
    return min(seq, rc)


def _all_canonical_6mers() -> list[str]:
    seen: set[str] = set()
    for combo in product("ACGT", repeat=6):
        seen.add(_canonical("".join(combo)))
    return sorted(seen)


def generate_priority_lists(
    n_landmarks: int = N_LANDMARKS,
    pinned_sites: tuple[str, ...] = DEFAULT_PRIORITY_SITES,
) -> list[list[str]]:
    """Build the per-landmark canonical-6-mer search-priority lists.

    The first `len(pinned_sites)` slots are anchored: slot `i` has
    `pinned_sites[i]` (canonicalized) at index 0 of its pile. Remaining
    slots get no pin — their pile starts with whatever round-robin
    distribution puts there, which spreads non-pinned slots away from any
    single MCS cluster. The remaining canonical 6-mers (minus the pinned
    set) are distributed round-robin across all piles after pinning, so no
    6-mer appears in more than one pile.
    """
    all_mers = _all_canonical_6mers()
    canonical_pinned = [_canonical(s) for s in pinned_sites[:n_landmarks]]
    pinned_set = set(canonical_pinned)
    pool = [m for m in all_mers if m not in pinned_set]

    piles: list[list[str]] = [[] for _ in range(n_landmarks)]
    for i, ps in enumerate(canonical_pinned):
        piles[i].append(ps)
    for i, mer in enumerate(pool):
        piles[i % n_landmarks].append(mer)
    return piles


# ---------- landmark search ----------

@dataclass
class LandmarkHit:
    """Forensic record of one landmark resolution."""
    slot: int
    site: Optional[str]       # canonical 6-mer matched, or None if all priorities exhausted
    position: Optional[int]   # position in the sequence
    upstream: Optional[str]   # 2 bases immediately upstream of `position`
    feature: Optional[int]    # 2-bit value derived from upstream[0] (the stored landmark)


def _find_all_occurrences(sequence: str, site: str) -> list[int]:
    """Return positions of all occurrences of `site` in either strand of `sequence`."""
    positions: set[int] = set()
    rc = str(Seq(site).reverse_complement())
    for pattern in {site, rc}:
        pos = 0
        while True:
            idx = sequence.find(pattern, pos)
            if idx == -1:
                break
            positions.add(idx)
            pos = idx + 1
    return sorted(positions)


def find_landmark(sequence: str, priority_list: list[str], slot: int) -> LandmarkHit:
    """Resolve one landmark slot against `sequence`. Return a LandmarkHit.

    Search algorithm:
      1. Walk the priority list top-to-bottom.
      2. For each candidate site, find all occurrences (forward + reverse complement).
      3. If any occurrences: pick one. Multiple hits → tiebreak alphabetically
         by the 10 bases immediately downstream of the site.
      4. Read 2 bases immediately upstream of the chosen position.
      5. Return the first of those 2 bases (2 bits) as the feature.

    If the chosen position is too close to the 5' end (< 2 bases upstream) we
    advance to the next candidate. If the entire priority list is exhausted,
    returns a LandmarkHit with site=None and feature=None.
    """
    seq_upper = sequence.upper()
    for site in priority_list:
        positions = _find_all_occurrences(seq_upper, site)
        if not positions:
            continue
        if len(positions) == 1:
            chosen = positions[0]
        else:
            chosen = min(
                positions,
                key=lambda p: seq_upper[p + 6:p + 16] if p + 16 <= len(seq_upper) else seq_upper[p + 6:],
            )
        if chosen < 2:
            # too close to the 5' end to read 2 upstream bases
            continue
        upstream = seq_upper[chosen - 2:chosen]
        if upstream[0] not in _BASE_TO_INT:
            continue  # upstream contains N or another non-canonical character
        return LandmarkHit(
            slot=slot,
            site=site,
            position=chosen,
            upstream=upstream,
            feature=_BASE_TO_INT[upstream[0]],
        )
    return LandmarkHit(slot=slot, site=None, position=None, upstream=None, feature=None)


def find_all_landmarks(
    sequence: str,
    priority_lists: Optional[list[list[str]]] = None,
) -> list[LandmarkHit]:
    """Resolve all 10 landmark slots in order. Returns a list of LandmarkHit."""
    if priority_lists is None:
        priority_lists = generate_priority_lists()
    return [find_landmark(sequence, priority_lists[i], slot=i) for i in range(len(priority_lists))]


def features_from_hits(hits: list[LandmarkHit], default: int = 0) -> list[int]:
    """Extract just the 2-bit feature values for storage in the barcode.
    Slots with no hit get `default` (the value is recoverable as an absence
    pattern at verify time when paired with the original landmark hits).
    """
    return [(h.feature if h.feature is not None else default) for h in hits]


# ---------- forensic comparison ----------

@dataclass
class LandmarkComparison:
    slot: int
    stored_feature: int
    found_hit: LandmarkHit
    site_present: bool          # was any priority site found at all?
    site_matches_stored: bool   # if we know the site at encode-time was X, did we find X?
    feature_matches: bool       # did the upstream-base value match what's stored?


def compare_landmarks(
    stored_features: list[int],
    sequence: str,
    priority_lists: Optional[list[list[str]]] = None,
    encoded_hits: Optional[list[LandmarkHit]] = None,
) -> list[LandmarkComparison]:
    """Compare stored landmark features against features computed from `sequence`.

    `encoded_hits` (if known — the encoder may have logged them at attest time)
    enables richer comparison. With just the stored 2-bit features in the
    barcode, "no hit at encode" and "feature was A" are indistinguishable;
    `encoded_hits` lets us tell them apart.
    """
    if priority_lists is None:
        priority_lists = generate_priority_lists()
    found = find_all_landmarks(sequence, priority_lists)
    out: list[LandmarkComparison] = []
    for slot, (stored, hit) in enumerate(zip(stored_features, found)):
        if encoded_hits is not None and slot < len(encoded_hits):
            enc_hit = encoded_hits[slot]
            site_matches = enc_hit.site == hit.site
            if enc_hit.site is None and hit.site is None:
                feature_matches = True            # both consistently absent
            elif enc_hit.site is None or hit.site is None:
                feature_matches = False           # one absent, one present
            else:
                feature_matches = hit.feature == stored
        else:
            site_matches = False
            feature_matches = hit.feature is not None and hit.feature == stored
        out.append(LandmarkComparison(
            slot=slot,
            stored_feature=stored,
            found_hit=hit,
            site_present=hit.site is not None,
            site_matches_stored=site_matches,
            feature_matches=feature_matches,
        ))
    return out
