# α-Lattice Design Randomization — Manim Animation

An animated explanation of the construction **and** randomization of α-designs
(alpha-lattice incomplete block designs), following Patterson & Williams (1976),
*A New Class of Resolvable Incomplete Block Designs*, Biometrika 63(1):83–92.

An introduction plus two parts, **14 scenes (S01–S14)** in playback order. Heat-map palette:
0 co-occurrences = blue, 1 = green, 2 = red, diagonal = grey; histograms count **unordered**
pairs (each pair once; totals = C(v,2)). ("FielDHub" is written **FieldHub** in
spoken/on-screen text so the text-to-speech pronounces it correctly.)
- **Intro (S01)** — the trade-off and what the design maximizes: small blocks lower the
  within-block variance σ², but cause confounding; the α-design maximizes the efficiency
  factor **E** by minimizing co-occurrences.
- **Part 1 (S02–S08)** — construction + randomization of the paper's Table 1 design:
  **v=20, r=3, s=4, k=5**. Here k>s, so some pairs are *forced* to co-occur twice → **α(0,1,2)**.
  Includes the co-occurrence heat-map (S07).
- **Part 2 (S09–S14)** — a **different design** (**v=30, r=3, s=6, k=5**) illustrating the
  improved **α(0,1)** lattice via the paper's Series III construction (same family as
  FieldHub). When **k ≤ s** no pair co-occurs more than once; the (abstract) improvement
  over Part 1's α(0,1,2) is the *lower co-occurrences*, which spread varieties more evenly
  (higher E) — not "more blocks of the same size," which would just be a bigger trial.
  Closes on the real **v=50** FieldHub trial.

Terminology: pairs **co-occur** (share a block); "concur" is avoided (it means *agree*).

## Layout

| Path | Tracked? | Purpose |
|------|----------|---------|
| `alpha_design.py`   | yes | Builds + verifies both designs; emits the two JSON files below |
| `design_data.json`  | yes | Part 1 data: v=20 α(0,1,2), provably matches Patterson & Williams Table 1 |
| `design_data_alpha01.json` | yes | Part 2 data: v=30 α(0,1) (Series III, k=5/s=6) + the real FieldHub v=50 trial census |
| `storyboard.md`     | yes | Scene-by-scene actions + voiceover narration |
| `alpha_lattice.py`  | yes | Manim scenes |
| `output/`           | **no** | Rendered MP4s (gitignored) |
| `media/`            | **no** | Manim cache/partial frames (gitignored) |
| `.venv/`            | **no** | Isolated Python 3.12 environment (gitignored) |

## Setup (run these yourself via `!` — they touch system/network)

```bash
brew install ffmpeg                          # manim hard-requires ffmpeg
cd alpha_lattice_animation
uv venv --python 3.12 .venv
uv pip install --python .venv -r requirements.txt
```

## Render

`build.sh` renders all 14 scenes and concatenates them into one narrated film at
`output/alpha_lattice_full.mp4` (macOS `say` narration is cached after the first run).

```bash
./build.sh         # fast 480p15 preview
./build.sh -qh     # final 1080p60
```

Render a single scene while iterating (per-scene file lands under `media/`, which is gitignored):

```bash
.venv/bin/manim -ql alpha_lattice.py S1_GeneratingArray
```

Optional: `brew install sox` silences the SoX warning and improves narration timing.
