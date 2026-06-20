"""Ground-truth data for the alpha-design animation, verified against the sources.

Part 1 - Patterson & Williams (1976), Table 1:  v=20, r=3, s=4, k=5.
  Here k=5 > s=4, so some pairs are FORCED to co-occur twice -> an alpha(0,1,2) design.

Part 2 - a small alpha(0,1) design built with the same Series III formula (Table 3)
  that FielDHub uses: v=12, r=3, s=4, k=3. Here k=3 <= s=4, so every pair co-occurs
  AT MOST once -> an alpha(0,1) design. This mirrors the user's real FielDHub trial
  (v=50, k=5, s=10, r=3), whose co-occurrence census is recorded in REAL_TRIAL.

The construction (paper section 3):
  A. Cyclic development: each column of the k x r generating array alpha spawns
     s columns by adding 1 (mod s) each step -> intermediate array alpha* (k x rs).
  B. Row offsets: add row_index * s to every entry -> variety labels 0..v-1.
  C. Columns are blocks; each group of s columns is a complete replication.

Pure stdlib (numpy not assumed present).
"""

import json
import os
from collections import Counter
from itertools import combinations


# --- generic construction --------------------------------------------------
def cyclic_development(alpha, s):
    """alpha (k x r) -> alpha_star (k x r*s) by adding 1 mod s, s times per column."""
    k, r = len(alpha), len(alpha[0])
    return [[(alpha[row][j] + t) % s for j in range(r) for t in range(s)]
            for row in range(k)]


def build_design(alpha, s):
    """Return (alpha_star, row_offsets, blocks, replications)."""
    k, r = len(alpha), len(alpha[0])
    star = cyclic_development(alpha, s)
    offsets = [row * s for row in range(k)]
    blocks = [[star[row][c] + offsets[row] for row in range(k)] for c in range(r * s)]
    reps = [blocks[j * s:(j + 1) * s] for j in range(r)]
    return star, offsets, blocks, reps


def cooccurrence_matrix(blocks, v):
    """Count, for each unordered pair of varieties, the number of shared blocks."""
    mat = [[0] * v for _ in range(v)]
    for block in blocks:
        for i, jj in combinations(sorted(block), 2):
            mat[i][jj] += 1
            mat[jj][i] += 1
    return mat


def offdiag_counts(mat, v):
    return Counter(mat[i][j] for i, j in combinations(range(v), 2))


def pairs_with(mat, v, value):
    return [[i, j] for i, j in combinations(range(v), 2) if mat[i][j] == value]


def blocks_with_pair(reps, a, b):
    """For pair (a,b): list of (rep_idx, block_idx, row_of_a, row_of_b) where both occur."""
    return [(ri, bi, blk.index(a), blk.index(b))
            for ri, rep in enumerate(reps) for bi, blk in enumerate(rep)
            if a in blk and b in blk]


def series_III(s, k):
    """Patterson & Williams Table 3, Series III generating array (r=3, s even, k<=s-1).
    Reduced form: column 1 zeros, column 2 the ramp 0..k-1, column 3 interleaves
    0, s/2, 1, s/2+1, 2, ...  (the construction FielDHub uses for r=3, even s)."""
    sh = s // 2
    col3 = [(i // 2) if i % 2 == 0 else sh + (i // 2) for i in range(k)]
    return [[0, i, col3[i]] for i in range(k)]


def emit(path, obj):
    with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), path), "w") as fh:
        json.dump(obj, fh, indent=2)


# ===========================================================================
# PART 1 - alpha(0,1,2): v=20, r=3, s=4, k=5  (Patterson & Williams Table 1)
# ===========================================================================
ALPHA20 = [
    [0, 0, 0],
    [0, 1, 2],
    [0, 2, 3],
    [0, 3, 1],
    [0, 3, 2],
]
STAR20, OFF20, BLOCKS20, REPS20 = build_design(ALPHA20, s=4)
CONC20 = cooccurrence_matrix(BLOCKS20, 20)

# assertions against the published paper
assert STAR20 == [
    [0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3],
    [0, 1, 2, 3, 1, 2, 3, 0, 2, 3, 0, 1],
    [0, 1, 2, 3, 2, 3, 0, 1, 3, 0, 1, 2],
    [0, 1, 2, 3, 3, 0, 1, 2, 1, 2, 3, 0],
    [0, 1, 2, 3, 3, 0, 1, 2, 2, 3, 0, 1],
], "alpha* does not match the paper"
assert REPS20 == [
    [[0, 4, 8, 12, 16], [1, 5, 9, 13, 17], [2, 6, 10, 14, 18], [3, 7, 11, 15, 19]],
    [[0, 5, 10, 15, 19], [1, 6, 11, 12, 16], [2, 7, 8, 13, 17], [3, 4, 9, 14, 18]],
    [[0, 6, 11, 13, 18], [1, 7, 8, 14, 19], [2, 4, 9, 15, 16], [3, 5, 10, 12, 17]],
], "blocks/replications do not match the paper"
assert CONC20[7][8] == 2, "varieties 7 and 8 should co-occur twice"
OFF20_COUNTS = offdiag_counts(CONC20, 20)
assert set(OFF20_COUNTS) == {0, 1, 2}, "part 1 should be alpha(0,1,2)"

# Tracked pair for the animation: co-occurs TWICE and, in BOTH shared blocks, the two
# varieties are non-adjacent (>=2 rows apart) -- so the student sees that "sharing a
# block" is about membership, not about being neighbouring plots.
PART1_TRACK = [4, 16]
_p1 = blocks_with_pair(REPS20, *PART1_TRACK)
assert len(_p1) == 2, "part 1 tracked pair must co-occur exactly twice"
assert all(abs(ra - rb) >= 2 for _, _, ra, rb in _p1), \
    "part 1 tracked pair must be non-adjacent within each shared block"

emit("design_data.json", {
    "params": {"v": 20, "r": 3, "s": 4, "k": 5},
    "alpha": ALPHA20,
    "alpha_star": STAR20,
    "row_offsets": OFF20,
    "replications": REPS20,
    "cooccurrence_offdiag_counts": {str(k): v for k, v in sorted(OFF20_COUNTS.items())},
    "cooccurrence_pairs_twice": pairs_with(CONC20, 20, 2),
    "cooccurrence_matrix": CONC20,
    "track_pair": PART1_TRACK,
    "design_class": "alpha(0, 1, 2)",
})


# ===========================================================================
# PART 2 - alpha(0,1) via Series III: v=30, r=3, s=6, k=5.
#   Same block size k=5 as Part 1, but now s=6 blocks per replication, so k <= s and no
#   pair is forced to co-occur twice. Uses Series III (r=3, s even) -- the SAME family of
#   construction as the FieldHub trial (k=5, s=10), so the two match.
# ===========================================================================
V_P2, R_P2, S_P2, K_P2 = 30, 3, 6, 5
ALPHA_P2 = series_III(s=S_P2, k=K_P2)     # [[0,0,0],[0,1,3],[0,2,1],[0,3,4],[0,4,2]]
STAR_P2, OFF_P2, BLOCKS_P2, REPS_P2 = build_design(ALPHA_P2, s=S_P2)
CONC_P2 = cooccurrence_matrix(BLOCKS_P2, V_P2)
OFF_P2_COUNTS = offdiag_counts(CONC_P2, V_P2)

# the whole point of part 2: k <= s, so NO pair co-occurs more than once
assert max(OFF_P2_COUNTS) == 1, "part 2 must be alpha(0,1): no pair co-occurs twice"
# tracked pair: co-occurs exactly once and is non-adjacent in its block (another plot
# sits between them) -- separating block co-occurrence from neighbouring-plot proximity.
TRACK_PAIR = None
for _i, _j in combinations(range(V_P2), 2):
    if CONC_P2[_i][_j] == 1:
        _bp = blocks_with_pair(REPS_P2, _i, _j)
        if abs(_bp[0][2] - _bp[0][3]) >= 2:
            TRACK_PAIR = [_i, _j]
            break
assert TRACK_PAIR is not None, "no non-adjacent single-co-occurrence pair found"

# --- the user's real FieldHub design (v=50, r=3, s=10, k=5), entries 1..50 ----
# Reproduced from FielDHub alpha_lattice(t=50,k=5,r=3,l=1,seed=1234); see rcbd_vs_alpha_lattice.qmd.
FH_BLOCKS = [
    [11, 23, 31, 37, 48], [1, 13, 17, 18, 24], [19, 33, 35, 41, 44], [16, 28, 30, 46, 49],
    [3, 10, 39, 40, 50], [21, 26, 29, 42, 43], [6, 8, 12, 14, 20], [2, 5, 25, 32, 36],
    [4, 27, 34, 38, 45], [7, 9, 15, 22, 47],
    [1, 2, 3, 19, 31], [12, 25, 37, 43, 49], [21, 30, 32, 33, 38], [23, 27, 40, 41, 42],
    [5, 10, 29, 35, 47], [13, 22, 34, 44, 46], [9, 14, 18, 45, 50], [6, 15, 16, 24, 36],
    [7, 11, 20, 26, 39], [4, 8, 17, 28, 48],
    [7, 17, 25, 33, 40], [4, 15, 19, 39, 43], [6, 10, 13, 37, 38], [1, 12, 27, 30, 47],
    [5, 18, 20, 23, 46], [2, 9, 26, 28, 41], [21, 36, 44, 48, 50], [3, 8, 22, 32, 42],
    [11, 24, 35, 45, 49], [14, 16, 29, 31, 34],
]
CONC_FH = cooccurrence_matrix([[e - 1 for e in blk] for blk in FH_BLOCKS], 50)  # 0-indexed
FH_COUNTS = offdiag_counts(CONC_FH, 50)
assert dict(FH_COUNTS) == {0: 925, 1: 300}, f"FieldHub design not alpha(0,1): {dict(FH_COUNTS)}"

emit("design_data_alpha01.json", {
    "params": {"v": V_P2, "r": R_P2, "s": S_P2, "k": K_P2},
    "alpha": ALPHA_P2,
    "alpha_star": STAR_P2,
    "row_offsets": OFF_P2,
    "replications": REPS_P2,
    "cooccurrence_offdiag_counts": {str(k): v for k, v in sorted(OFF_P2_COUNTS.items())} | {"2": 0},
    "cooccurrence_matrix": CONC_P2,
    "track_pair": TRACK_PAIR,
    "design_class": "alpha(0, 1)",
    # the user's real design, generated by FielDHub (read aloud as "FieldHub"):
    "real_trial": {
        "software": "FieldHub", "series": "III",
        "params": {"v": 50, "r": 3, "s": 10, "k": 5},
        "cooccurrence_offdiag_counts": {str(k): v for k, v in sorted(FH_COUNTS.items())} | {"2": 0},
        "cooccurrence_matrix": CONC_FH,
        "n_pairs": 1225,
        "design_class": "alpha(0, 1)",
    },
})


# --- report ----------------------------------------------------------------
if __name__ == "__main__":
    print("All assertions passed.\n")
    print("PART 1  alpha(0,1,2)  v=20 r=3 s=4 k=5")
    print(f"  off-diagonal co-occurrence counts: {dict(sorted(OFF20_COUNTS.items()))}")
    print(f"  pairs co-occurring twice: {len(pairs_with(CONC20, 20, 2))}")
    print("\nPART 2  alpha(0,1)    v=30 r=3 s=6 k=5   (Series III; matches FieldHub family)")
    print(f"  generating array: {ALPHA_P2}")
    for i, rep in enumerate(REPS_P2, 1):
        print(f"  Rep {i}: " + "  ".join(str(b) for b in rep))
    print(f"  off-diagonal co-occurrence counts: {dict(sorted(OFF_P2_COUNTS.items()))}")
    print(f"  tracked pair (co-occurs once): {TRACK_PAIR}")
    print("\nWrote design_data.json and design_data_alpha01.json")
