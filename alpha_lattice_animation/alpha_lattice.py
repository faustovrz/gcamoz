"""Manim animation: construction and randomization of an alpha-lattice design.

Follows Patterson & Williams (1976), Table 1 (v=20, r=3, s=4, k=5).
All numbers are read from design_data.json (produced by alpha_design.py) so the
animation can never drift from the verified ground truth.

Narration uses the macOS `say` voice via the local say_service. Scenes are S01..S14,
in playback order; build and concatenate them with ./build.sh (see README).
"""

import json
import os
import random
import sys

from manim import (
    BLUE,
    DOWN,
    GOLD,
    GREEN,
    GREY,
    LEFT,
    MAROON,
    ORIGIN,
    RED,
    RIGHT,
    UP,
    WHITE,
    Arrow,
    Brace,
    Create,
    FadeIn,
    FadeOut,
    Indicate,
    Integer,
    MathTex,
    Rectangle,
    ReplacementTransform,
    Scene,
    Square,
    SurroundingRectangle,
    Text,
    Transform,
    VGroup,
    Write,
)
from manim_voiceover import VoiceoverScene

# --- load verified design data --------------------------------------------
_HERE = os.path.dirname(os.path.abspath(__file__))
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)  # so `say_service` resolves when manim loads this file
from say_service import SayService  # noqa: E402  (local macOS `say` TTS service)
with open(os.path.join(_HERE, "design_data.json")) as _fh:
    D = json.load(_fh)
with open(os.path.join(_HERE, "design_data_alpha01.json")) as _fh:
    D2 = json.load(_fh)  # Part 2: the α(0,1) design + the real FielDHub trial summary

V, R, S, K = (D["params"][x] for x in ("v", "r", "s", "k"))
ALPHA = D["alpha"]
ALPHA_STAR = D["alpha_star"]
ROW_OFFSETS = D["row_offsets"]
REPS = D["replications"]

# --- shared visual config (define colors once) -----------------------------
C_SEED = BLUE          # generating-array entries
C_OFFSET = GOLD        # row offset amounts
C_HILITE = MAROON      # the varieties 7 & 8 co-occurrence highlight
REP_COLORS = [BLUE, GREEN, GOLD]   # one per replication
CELL = 0.62            # cell side length
FS = 26                # number font size

# co-occurrence heat-map palette: 0 -> blue, 1 -> green, 2 -> red; diagonal (self) -> grey
HEAT = {0: BLUE, 1: GREEN, 2: RED}
HEAT_DIAG = GREY


class NumberGrid(VGroup):
    """A grid of square cells holding integers. Access cells via self.cell(r, c)."""

    def __init__(self, values, color=WHITE, cell=CELL, fs=FS, **kwargs):
        super().__init__(**kwargs)
        self.rows = len(values)
        self.cols = len(values[0])
        self._cells = [[None] * self.cols for _ in range(self.rows)]
        for r in range(self.rows):
            for c in range(self.cols):
                sq = Square(side_length=cell, color=color, stroke_width=2)
                sq.move_to([c * cell, -r * cell, 0])
                num = Integer(values[r][c], font_size=fs).move_to(sq.get_center())
                grp = VGroup(sq, num)
                self._cells[r][c] = grp
                self.add(grp)
        self.center()

    def cell(self, r, c):
        return self._cells[r][c]

    def square(self, r, c):
        return self._cells[r][c][0]

    def number(self, r, c):
        return self._cells[r][c][1]

    def column(self, c):
        return VGroup(*(self._cells[r][c] for r in range(self.rows)))

    def row(self, r):
        return VGroup(*(self._cells[r][c] for c in range(self.cols)))


def make_heatmap(matrix, cell):
    """A v x v grid of colored squares: HEAT[value], diagonal grey. No strokes (fast)."""
    v = len(matrix)
    grid = VGroup()
    for i in range(v):
        for j in range(v):
            color = HEAT_DIAG if i == j else HEAT.get(matrix[i][j], RED)
            sq = Square(side_length=cell, stroke_width=0)
            sq.set_fill(color, opacity=1.0)
            sq.move_to([j * cell, -i * cell, 0])
            grid.add(sq)
    grid.center()
    return grid


def make_legend(font_size=22):
    items = [("0  never", HEAT[0]), ("1  once", HEAT[1]), ("2  twice", HEAT[2])]
    row = VGroup()
    for label, color in items:
        sw = Square(0.28, stroke_width=0).set_fill(color, 1)
        row.add(VGroup(sw, Text(label, font_size=font_size).next_to(sw, RIGHT, buff=0.12))
                .arrange(RIGHT, buff=0.12))
    return row.arrange(RIGHT, buff=0.5)


def make_histogram(counts, height=2.6, bar_w=1.0, gap=0.7):
    """Labeled bar chart from {0:n0, 1:n1, 2:n2}; bar color = HEAT[k]; count above, label below."""
    keys = sorted(int(k) for k in counts)
    vals = [int(counts[str(k)]) for k in keys]
    maxc = max(vals) or 1
    chart = VGroup()
    for idx, k in enumerate(keys):
        h = max(height * vals[idx] / maxc, 0.04)
        bar = Rectangle(width=bar_w, height=h, stroke_color=WHITE, stroke_width=1.5)
        bar.set_fill(HEAT[k], opacity=1.0)
        bar.move_to([idx * (bar_w + gap), h / 2, 0])
        cnt = Text(str(vals[idx]), font_size=26, weight="BOLD").next_to(bar, UP, buff=0.12)
        xlab = Text(f"{k}x", font_size=24, color=HEAT[k]).next_to(bar, DOWN, buff=0.18)
        chart.add(VGroup(bar, cnt, xlab))
    return chart.center()


def alpha_label(text, font_size=28):
    """White 'alpha(...) design' label with the digits coloured to match the heat-map
    palette (0 blue, 1 green, 2 red)."""
    return Text(text, font_size=font_size, color=WHITE,
                t2c={"0": HEAT[0], "1": HEAT[1], "2": HEAT[2]})


def make_design_preview(r, s, k, cell=0.34):
    """Empty, colour-by-rep preview of the final design (no numbers): r reps, each s blocks
    of k plots. Plots within a block touch (no gap); gaps between blocks; larger between reps."""
    reps = VGroup()
    for ri in range(r):
        blocks = VGroup()
        for _ in range(s):
            col = VGroup(*[
                Square(cell, stroke_width=2).set_stroke(REP_COLORS[ri % len(REP_COLORS)])
                .set_fill(REP_COLORS[ri % len(REP_COLORS)], 0.30)
                for _ in range(k)
            ]).arrange(DOWN, buff=0)          # no space between plots in a block
            blocks.add(col)
        blocks.arrange(RIGHT, buff=0.14)       # gap between blocks
        reps.add(blocks)
    reps.arrange(RIGHT, buff=0.6)              # larger gap between reps
    return reps


def voice(scene):
    """Attach the macOS `say` voiceover service (natural, local) to a scene."""
    scene.set_speech_service(SayService(
        voice="Evan (Enhanced)", rate=155, sentence_pause=550, comma_pause=180,
        dash_pause=320, trailing_pause=750, letter_pause=250))


# ===========================================================================
class S01_Optimize(VoiceoverScene):
    """Introduction: the trade-off, and what the alpha-design actually maximizes (E)."""

    def construct(self):
        voice(self)
        title = Text("What an alpha-design optimizes", font_size=40,
                     weight="BOLD").to_edge(UP)
        # Opening title card -- exactly the S02 card (centered name + citation), shown from the
        # start so "the paper" has a referent, until the box illustration appears.
        card = VGroup(
            Text("Alpha-Lattice Designs", font_size=48, weight="BOLD"),
            Text("Patterson & Williams (1976)", font_size=24),
        )
        card[1].next_to(card[0], DOWN)
        self.add(card)
        self.wait(2.5)   # silent title card (no voiceover) before the content begins

        # a field strip with a fertility trend (opacity ramp = soil gradient)
        field = VGroup()
        for i in range(10):
            sq = Square(side_length=0.7)
            sq.set_fill(GREEN, opacity=0.12 + 0.085 * i)
            sq.set_stroke(WHITE, 1)
            field.add(sq)
        field.arrange(RIGHT, buff=0.06).move_to(ORIGIN)
        grad_label = Text("a field with a fertility trend across it",
                          font_size=24).next_to(field, DOWN, buff=0.4)

        with self.voiceover(
            text="Why split a trial into small incomplete blocks at all? It comes down to "
            "one trade-off. Picture a field with a fertility trend running across it."
        ):
            self.play(FadeOut(card), FadeIn(title), FadeIn(field), FadeIn(grad_label))

        big = SurroundingRectangle(field, color=GOLD, buff=0.12)
        big_lbl = Text("one big block: large within-block variance",
                       font_size=24, color=GOLD).next_to(title, DOWN, buff=0.2)
        with self.voiceover(
            text="If a whole replication is one big block, plots at opposite ends sit on "
            "very different soil, so the within-block error variance is large, and real "
            "variety differences are hard to see."
        ):
            self.play(Create(big), FadeIn(big_lbl))

        smalls = VGroup(*(
            SurroundingRectangle(VGroup(field[2 * j], field[2 * j + 1]),
                                 color=BLUE, buff=0.08)
            for j in range(5)))
        small_lbl = Text("small blocks: nearly uniform soil, small within-block variance",
                         font_size=24, color=BLUE).next_to(title, DOWN, buff=0.2)
        with self.voiceover(
            text="Break it into small blocks and each block covers nearly uniform soil. "
            "The within-block variance drops, and comparisons within a block get sharper."
        ):
            self.play(ReplacementTransform(big, smalls),
                      ReplacementTransform(big_lbl, small_lbl))
        self.play(FadeOut(VGroup(field, grad_label, smalls, small_lbl)))

        formula = MathTex(
            r"\text{precision of a comparison}\;\propto\;\frac{1}{\sigma^2}\,\times\,E",
            font_size=44).move_to(ORIGIN)
        sigma_line = MathTex(r"\sigma^2 = \text{within-block error variance}",
                             font_size=28, color=WHITE)
        eff_line = Text("efficiency  E = fraction of precision kept,  0 < E ≤ 1",
                        font_size=26, color=C_HILITE)
        min_line = Text("by minimizing co-occurrences", font_size=26)
        expl = VGroup(sigma_line, eff_line, min_line).arrange(DOWN, buff=0.3)
        expl.next_to(formula, DOWN, buff=0.7)

        with self.voiceover(
            text="But there is a catch: varieties in different blocks are not compared "
            "directly, so some information is lost. Precision depends on two things -- the "
            "within-block variance sigma squared, and the efficiency E."
        ):
            self.play(Write(formula))
            self.play(FadeIn(sigma_line, shift=UP * 0.2))

        with self.voiceover(
            text="The block size sets the within-block variance, and that is the "
            "experimenter's choice."
        ):
            pass

        # red efficiency line appears as the narration starts the definition
        with self.voiceover(
            text="The efficiency factor E, between zero and one, is the fraction of "
            "comparison precision you keep despite the blocking."
        ):
            self.play(FadeIn(eff_line, shift=UP * 0.2))

        with self.voiceover(text="The design maximizes that efficiency in just one way."):
            pass

        # "by minimizing co-occurrences" line appears as the narration reaches that point
        with self.voiceover(
            text="By minimizing co-occurrences -- keeping any two varieties from sharing a "
            "block more often than they must. That one goal drives everything that follows."
        ):
            self.play(FadeIn(min_line, shift=UP * 0.2))
        self.play(FadeOut(VGroup(title, formula, expl)))


# ===========================================================================
class S02_Problem(VoiceoverScene):
    def construct(self):
        voice(self)
        title = Text("Alpha-Lattice Designs", font_size=44, weight="BOLD").to_edge(UP)
        sub = Text("Patterson & Williams (1976)", font_size=24).next_to(title, DOWN, buff=0.12)

        # Preview of the final design: colour-by-rep, split block-cards, no numbers.
        grid = make_design_preview(R, S, K)
        grid.move_to(ORIGIN).shift(DOWN * 0.25)
        k_brace = Brace(grid[0][0], LEFT)
        k_lbl = k_brace.get_text(f"k = {K}")
        s_brace = Brace(grid[0], UP)
        s_lbl = s_brace.get_text(f"s = {S}")
        r_brace = Brace(grid, DOWN)
        r_lbl = r_brace.get_text(f"r = {R}")
        v_note = MathTex(rf"v = k \times s = {V}", font_size=34).next_to(r_lbl, DOWN, buff=0.3)

        with self.voiceover(
            text="Suppose we want to compare twenty varieties in a field trial, "
            "with three full replications."
        ):
            self.play(Write(title))
            self.play(FadeIn(sub, shift=UP * 0.2))

        with self.voiceover(
            text="A whole replication in one long strip is noisy, because distant plots "
            "differ in soil and moisture. So we break each replication into smaller, more "
            "uniform incomplete blocks: four blocks of five plots each."
        ):
            self.play(FadeIn(grid))
            self.play(FadeIn(k_brace), FadeIn(k_lbl), FadeIn(s_brace), FadeIn(s_lbl))

        with self.voiceover(
            text="Three replications in all, and twenty varieties, since twenty is five "
            "times four. The challenge is choosing which varieties share a block, so that "
            "comparisons stay balanced. Patterson and Williams gave a simple recipe."
        ):
            self.play(FadeIn(r_brace), FadeIn(r_lbl), FadeIn(v_note))
            self.play(Indicate(grid, scale_factor=1.03))
        self.play(FadeOut(VGroup(title, sub, grid, k_brace, k_lbl, s_brace, s_lbl,
                                 r_brace, r_lbl, v_note)))


# ===========================================================================
class S03_GeneratingArray(VoiceoverScene):
    def construct(self):
        voice(self)
        grid = NumberGrid(ALPHA, color=C_SEED, cell=0.5, fs=22)  # cell 0.5 matches S04's alpha*
        label = Text("Generating array  (5 rows x 3 columns)", font_size=28)
        label.next_to(grid, UP, buff=0.6)
        chip = Text("k = 5      s = 4      r = 3", font_size=26).to_edge(DOWN)

        with self.voiceover(
            text="The recipe starts with a generating array: five rows for the five plots "
            "in a block, three columns for the three replications."
        ):
            self.play(Write(label))
            self.play(Create(grid), FadeIn(chip))

        with self.voiceover(
            text="Every entry is a number from zero to three, a residue modulo four, "
            "where four is the number of blocks in a replication."
        ):
            mod_note = VGroup(
                Text("residue modulo 4", font_size=26),
                Text("(4 = s, the number of blocks/rep)", font_size=22, color=C_OFFSET),
            ).arrange(DOWN, buff=0.15).next_to(grid, RIGHT, buff=0.9)
            self.play(FadeIn(mod_note), Indicate(grid, scale_factor=1.05))

        with self.voiceover(
            text="This array is in reduced form: the whole first row and first column "
            "are zero. These few numbers are the entire seed of the design."
        ):
            row0 = SurroundingRectangle(grid.row(0), color=C_OFFSET, buff=0.05)
            col0 = SurroundingRectangle(grid.column(0), color=C_OFFSET, buff=0.05)
            self.play(Create(row0), Create(col0))
        self.play(FadeOut(VGroup(row0, col0)))

        with self.voiceover(
            text="And these numbers are not arbitrary. This little array is precomputed -- "
            "looked up from tables, given by a formula, or found by search -- and chosen so "
            "that the finished design maximizes the efficiency factor E. Picking the array "
            "is the only clever step; everything after it is mechanical."
        ):
            note = Text("precomputed to maximize E", font_size=26, color=C_HILITE)
            note.next_to(grid, DOWN, buff=0.5)
            self.play(Write(note), Indicate(grid, scale_factor=1.05))
        self.play(FadeOut(VGroup(label, note, chip, mod_note)))
        self.grid = grid  # for visual continuity if rendered together


# ===========================================================================
class S04_CyclicDevelopment(VoiceoverScene):
    def construct(self):
        voice(self)
        # Open on S03's three generating-array columns (contiguous, centered), then split
        # them apart into the alpha* layout. cell=0.5 matches S03 for a seamless cut.
        star = NumberGrid(ALPHA_STAR, color=C_SEED, cell=0.5, fs=22)
        seed = [0, 4, 8]                       # alpha* columns that equal the generating array
        target_x = {c: star.column(c).get_x() for c in seed}
        for i, c in enumerate(seed):           # move seeds to contiguous, centered positions
            star.column(c).shift(RIGHT * ((i - 1) * 0.5 - target_x[c]))
        seed_group = VGroup(*(star.column(c) for c in seed))
        self.add(seed_group)
        label = Text("Cyclic development  (+1 mod 4, wrapping 3 -> 0)", font_size=26)
        label.next_to(seed_group, UP, buff=0.6)

        with self.voiceover(
            text="Now we develop each column cyclically. First, the three columns of the "
            "generating array slide apart, to make room."
        ):
            self.play(Write(label))
            self.play(*[
                star.column(c).animate.shift(RIGHT * (target_x[c] - star.column(c).get_x()))
                for c in seed
            ])

        with self.voiceover(
            text="Take each column and make three more copies, adding one to every number "
            "each time -- and whenever a three would become four, it wraps back around to zero."
        ):
            self.play(*[Create(star.column(c)) for c in (1, 2, 3)])

        with self.voiceover(
            text="Each original column blossoms into a group of four. The three columns of "
            "the seed become twelve columns in this intermediate array."
        ):
            self.play(*[Create(star.column(c)) for c in (5, 6, 7, 9, 10, 11)])
        self.play(FadeOut(label))
        # Leave the grid centered (no upward move) so it merges seamlessly with S05.


# ===========================================================================
class S05_RowOffsets(VoiceoverScene):
    def construct(self):
        voice(self)
        star = NumberGrid(ALPHA_STAR, color=C_SEED, cell=0.5, fs=22)
        offsets = VGroup(
            *[Integer(o, font_size=24, color=C_OFFSET) for o in ROW_OFFSETS]
        )
        for r in range(K):
            offsets[r].next_to(star.row(r), LEFT, buff=0.5)
        plus = VGroup(*[Text("+", font_size=24, color=C_OFFSET).next_to(offsets[r], RIGHT, buff=0.1) for r in range(K)])

        self.add(star)
        with self.voiceover(
            text="So far every number is small, between zero and three. To turn them into "
            "variety labels from zero to nineteen, we shift each row by a fixed amount."
        ):
            self.play(FadeIn(offsets, shift=RIGHT * 0.2), FadeIn(plus))

        with self.voiceover(
            text="The first row stays. We add four to the second row, eight to the third, "
            "twelve to the fourth, and sixteen to the fifth."
        ):
            anims = []
            for r in range(K):
                for c in range(R * S):
                    new_val = ALPHA_STAR[r][c] + ROW_OFFSETS[r]
                    target = Integer(new_val, font_size=22).move_to(star.number(r, c))
                    anims.append(Transform(star.number(r, c), target))
            self.play(*anims, run_time=2.0)

        with self.voiceover(
            text="Now all twenty varieties appear, and each one shows up exactly once in "
            "every group of four columns."
        ):
            self.play(Indicate(star, scale_factor=1.03))
        self.play(FadeOut(VGroup(offsets, plus)))


# ===========================================================================
class S06_Blocks(VoiceoverScene):
    def construct(self):
        voice(self)
        # Open on the exact blue variety matrix that S05 ended on (seamless cut).
        bm = _block_matrix(D)                       # 5 x 12 varieties; columns are blocks
        grid = NumberGrid(bm, color=C_SEED, cell=0.5, fs=22)
        self.add(grid)
        groups = [list(range(j * S, (j + 1) * S)) for j in range(R)]  # columns per rep

        with self.voiceover(
            text="These twelve columns are the blocks of the design. First, watch them "
            "separate at the replication boundaries."
        ):
            shifts = [LEFT * 0.9, ORIGIN, RIGHT * 0.9]
            self.play(*[
                VGroup(*(grid.column(c) for c in cols)).animate.shift(shifts[gi])
                for gi, cols in enumerate(groups)
            ])

        # Colour and label each replication in turn, each synced to naming it.
        ordinals = ["one", "two", "three"]
        labels = VGroup()
        for gi, cols in enumerate(groups):
            grp = VGroup(*(grid.column(c) for c in cols))
            lbl = Text(f"Rep {gi + 1}", font_size=26, color=REP_COLORS[gi]).next_to(grp, UP, buff=0.3)
            with self.voiceover(text=f"Replication {ordinals[gi]}."):
                self.play(
                    *[grid.square(r, c).animate.set_stroke(REP_COLORS[gi], width=3)
                      for c in cols for r in range(K)],
                    FadeIn(lbl, shift=DOWN * 0.2),
                )
            labels.add(lbl)

        with self.voiceover(
            text="Each group of four columns is one complete replication, containing all "
            "twenty varieties. So the design is resolvable."
        ):
            self.play(Indicate(labels, scale_factor=1.05))

        # Split each replication's four columns into separated block-cards. The 0.18 gap
        # matches the within-rep block spacing used in the randomization section (S08).
        with self.voiceover(
            text="And each column is a block -- a set of five varieties grown together."
        ):
            self.play(*[
                grid.column(c).animate.shift(RIGHT * ((c % S) - (S - 1) / 2) * 0.18)
                for c in range(R * S)
            ])

        # tracked pair (4 & 16): co-occurs twice, non-adjacent within each shared column
        PA, PB = D["track_pair"]
        rects = VGroup()
        for c in range(R * S):
            col_vals = [bm[r][c] for r in range(K)]
            if PA in col_vals and PB in col_vals:
                for val in (PA, PB):
                    rects.add(SurroundingRectangle(grid.cell(col_vals.index(val), c),
                                                   color=C_HILITE, buff=0.03, stroke_width=4))
        with self.voiceover(
            text="Look at varieties four and sixteen. They share a block twice -- once in "
            "replication one and once in replication three. And notice they are not "
            "neighbours: other plots sit between them in the column. Co-occurring means "
            "sharing a block, not being side by side."
        ):
            self.play(Create(rects))

        with self.voiceover(
            text="Most pairs meet zero or one times; a few meet twice. That tight "
            "co-occurrence structure makes the design efficient. This is an alpha zero, "
            "one, two design."
        ):
            cls = Text("α(0, 1, 2) design", font_size=30, color=C_HILITE).to_edge(DOWN)
            self.play(Write(cls))
        self.play(FadeOut(VGroup(grid, labels, cls, rects)))


# ===========================================================================
class S07_CoocPart1(VoiceoverScene):
    """Co-occurrence heat-map + histogram for the v=20 α(0,1,2) design."""

    def construct(self):
        voice(self)
        title = Text("Co-occurrence structure  (20 varieties)", font_size=32).to_edge(UP)
        self.add(title)

        heat = make_heatmap(D["cooccurrence_matrix"], cell=0.26)
        heat.to_edge(LEFT, buff=1.0).shift(DOWN * 0.3)
        legend = make_legend().next_to(heat, DOWN, buff=0.35)
        with self.voiceover(
            text="Here is the whole co-occurrence structure as a heat map. Each cell is a "
            "pair of varieties: blue if they never share a block, green if they share one, "
            "red if they share two. The grey diagonal is a variety with itself."
        ):
            self.play(FadeIn(heat), FadeIn(legend))

        hist = make_histogram(D["cooccurrence_offdiag_counts"]).scale(0.9)
        hist.to_edge(RIGHT, buff=1.0).shift(DOWN * 0.3)
        with self.voiceover(
            text="Counting them: eighty-two pairs never meet, ninety-six meet once, and "
            "twelve meet twice."
        ):
            self.play(FadeIn(hist, shift=UP * 0.2))

        acls = alpha_label("α(0, 1, 2) design").next_to(hist, DOWN, buff=0.45)
        with self.voiceover(
            text="Those twelve red cells are what make this an alpha, zero one two, design."
        ):
            self.play(Write(acls))
        self.wait(0.3)
        self.play(FadeOut(VGroup(title, heat, legend, hist, acls)))


# ===========================================================================
class S08_Randomization(VoiceoverScene):
    def construct(self):
        voice(self)
        # Independent randomization at every level (seeded so renders are reproducible):
        # each replication, each block within a rep, and each block's plots are shuffled
        # on their own -- the plot shuffle is NOT shared across blocks.
        rng = random.Random(20260619)
        rep_order = rng.sample(range(R), R)                          # 1. order of reps
        block_orders = [rng.sample(range(S), S) for _ in range(R)]   # 2. blocks within each rep
        plot_orders = [[rng.sample(range(K), K) for _ in range(S)]   # 3. plots within each block
                       for _ in range(R)]
        relabel = {old: (old * 3 + 7) % V for old in range(V)}       # 4. coprime relabel
        TA, TB = D["track_pair"]              # tracked pair: co-occurs twice, non-adjacent

        def make_rep_group(reps):
            rep_groups = VGroup()
            for ri, rep in enumerate(reps):
                blocks = VGroup()
                for blk in rep:
                    col = NumberGrid([[x] for x in blk], color=REP_COLORS[ri % R],
                                     cell=0.5, fs=22)
                    blocks.add(col)
                blocks.arrange(RIGHT, buff=0.18)
                rep_groups.add(blocks)
            rep_groups.arrange(RIGHT, buff=0.7).move_to(ORIGIN)  # cell 0.5, matches S06
            return rep_groups

        def pair_boxes(group, data, a, b):
            """Maroon boxes around every cell holding a or b in a block that has both."""
            boxes = VGroup()
            for ri, rep in enumerate(data):
                for bi, blk in enumerate(rep):
                    if a in blk and b in blk:
                        for val in (a, b):
                            row = blk.index(val)
                            boxes.add(SurroundingRectangle(
                                group[ri][bi].cell(row, 0),
                                color=C_HILITE, buff=0.05, stroke_width=4))
            return boxes

        def make_badge(a, b):
            return VGroup(
                Text(f"varieties {a} & {b} share a block:", font_size=22),
                Text("2x", font_size=26, color=C_HILITE, weight="BOLD"),
            ).arrange(RIGHT, buff=0.2).to_edge(DOWN)

        data = [[list(blk) for blk in rep] for rep in REPS]
        group = make_rep_group(data)
        boxes = pair_boxes(group, data, TA, TB)
        badge = make_badge(TA, TB)
        title = Text("Randomize before planting", font_size=30).to_edge(UP)
        self.add(title, group)

        with self.voiceover(
            text="The plan is systematic, so before planting we randomize it. But watch "
            "varieties four and sixteen, which share a block twice -- and note they are "
            "never side by side; other plots sit between them. The maroon boxes will follow "
            "them through every shuffle."
        ):
            self.play(Create(boxes), FadeIn(badge))

        def step(new_data, a, b, narration):
            nonlocal group, boxes
            new_group = make_rep_group(new_data)
            new_boxes = pair_boxes(new_group, new_data, a, b)
            with self.voiceover(text=narration):
                self.play(ReplacementTransform(group, new_group),
                          Transform(boxes, new_boxes),
                          Indicate(badge[1]))
            group = new_group

        # 1. permute replications -- slide WHOLE reps to new slots (numbers kept intact),
        #    then recolour them by their new position.
        d1 = [data[i] for i in rep_order]
        slot_x = [group[i].get_x() for i in range(R)]
        dest_x = [slot_x[rep_order.index(i)] for i in range(R)]   # old rep i -> its new slot
        moved_boxes = boxes.copy()
        for nb in moved_boxes:                                    # boxes ride along with reps
            j = min(range(R), key=lambda q: abs(nb.get_x() - slot_x[q]))
            nb.shift(RIGHT * (dest_x[j] - slot_x[j]))
        with self.voiceover(
            text="First we shuffle the order of the three replications. Whole replications "
            "switch places -- the numbers inside each one stay exactly as they were."
        ):
            self.play(*[group[i].animate.set_x(dest_x[i]) for i in range(R)],
                      Transform(boxes, moved_boxes), Indicate(badge[1]))
        # smooth recolour to the new positions -- no narration, just a colour transition
        self.play(*[
            group[i][bi].square(rr, 0).animate.set_stroke(REP_COLORS[rep_order.index(i)])
            for i in range(R) for bi in range(S) for rr in range(K)
        ], run_time=1.2)
        # swap in canonical-order mobjects (identical on screen) so later steps line up
        self.remove(group, boxes)
        group = make_rep_group(d1)
        boxes = pair_boxes(group, d1, TA, TB)
        self.add(group, boxes)

        # 2. permute blocks within each replication (independent order per rep)
        d2 = [[rep[bk] for bk in block_orders[i]] for i, rep in enumerate(d1)]
        step(d2, TA, TB,
             "Then, within each replication, we shuffle the order of the blocks. Two "
             "boxes still, two shared blocks still.")

        # 3. permute plots within each block (independent order per block)
        d3 = [[[blk[p] for p in plot_orders[i][j]] for j, blk in enumerate(rep)]
              for i, rep in enumerate(d2)]
        step(d3, TA, TB,
             "Within each block, we shuffle the plot positions. Four and sixteen slide "
             "around inside their blocks -- the gap between them even changes -- but they "
             "stay in the same block. Proximity is not the same as sharing a block.")

        # 4. relabel varieties (boxes stay put: same blocks, new names)
        d4 = [[[relabel[x] for x in blk] for blk in rep] for rep in d3]
        na, nb = relabel[TA], relabel[TB]
        new_group = make_rep_group(d4)
        new_boxes = pair_boxes(new_group, d4, na, nb)
        new_badge = make_badge(na, nb)
        with self.voiceover(
            text="Finally we relabel every variety. Four becomes nineteen, and sixteen "
            "becomes fifteen. The boxes do not move, because it is the very same two "
            "blocks. Only the names changed."
        ):
            self.play(ReplacementTransform(group, new_group),
                      Transform(boxes, new_boxes))
            self.play(Transform(badge, new_badge))
        group = new_group

        with self.voiceover(
            text="The co-occurrence count never changed: still exactly two. Randomization "
            "preserves the design's co-occurrence structure, while spreading any field "
            "trend randomly across varieties. A fair, efficient, ready-to-plant "
            "alpha-lattice layout."
        ):
            self.play(Indicate(boxes), Indicate(badge))
            done = Text("Construct once, randomize per site.", font_size=28)
            done.next_to(badge, UP, buff=0.3)
            self.play(Write(done))
        self.wait(0.5)


# ===========================================================================
# PART 2 -- when k ≤ s, the α(0,1) design: the paper's dedicated construction
#           (Series generators) caps every pair at one co-occurrence.
# ===========================================================================
def _block_matrix(d):
    """Developed + row-offset matrix (k rows x r*s cols); each column is a block."""
    star, off = d["alpha_star"], d["row_offsets"]
    cols = len(star[0])
    return [[star[i][c] + off[i] for c in range(cols)] for i in range(len(star))]


# Part 2 grid geometry (used by S11 to reproduce the design layout as numbered block-cards).
P2_CELL, P2_FS, P2_REPGAP, P2_COLGAP = 0.5, 20, 0.9, 0.14


def _split_shift(c, r, s, repgap, colgap):
    """Total x-shift moving column c of a centered r*s grid into its split block-card
    position: replications separated by `repgap`, columns within a rep gapped by `colgap`
    (symmetric, so each rep expands around its own centre)."""
    return (c // s - (r - 1) / 2) * repgap + ((c % s) - (s - 1) / 2) * colgap


class S09_Alpha01Design(VoiceoverScene):
    """Part 2 design overview, Scene-2 style: colour-by-rep block-cards, no numbers, braces."""

    def construct(self):
        voice(self)
        p = D2["params"]
        r, s, k, v = p["r"], p["s"], p["k"], p["v"]
        title = Text("Part 2  -  the α(0,1) design", font_size=34,
                     weight="BOLD").to_edge(UP)
        subtitle = Text("when k ≤ s, no pair co-occurs more than once   (here  k = 5,  s = 6)",
                        font_size=26, color=C_HILITE).next_to(title, DOWN, buff=0.15)
        self.add(title, subtitle)

        grid = make_design_preview(r, s, k, cell=0.30)
        grid.move_to(ORIGIN).shift(DOWN * 0.35)
        k_brace = Brace(grid[0][0], LEFT)
        k_lbl = k_brace.get_text(f"k = {k}")
        s_brace = Brace(grid[0], UP)
        s_lbl = s_brace.get_text(f"s = {s}")
        r_brace = Brace(grid, DOWN)
        r_lbl = r_brace.get_text(f"r = {r}")
        v_note = MathTex(rf"v = k \times s = {v}", font_size=32).next_to(r_lbl, DOWN, buff=0.28)

        with self.voiceover(
            text="Part two is a different design -- a separate trial, with thirty varieties -- "
            "illustrating the improved lattice you get when the block size k is at most the "
            "number of blocks s. Here, blocks of five plots and six blocks per replication: k "
            "equals five, s equals six."
        ):
            self.play(FadeIn(grid))
            self.play(FadeIn(k_brace), FadeIn(k_lbl), FadeIn(s_brace), FadeIn(s_lbl))
        with self.voiceover(
            text="Three replications -- thirty varieties in all, since thirty is five times "
            "six. Because k is at most s, no pair is ever forced to share a block more than "
            "once: the improvement of an alpha, zero one, design over an alpha zero one two."
        ):
            self.play(FadeIn(r_brace), FadeIn(r_lbl), FadeIn(v_note))
            self.play(Indicate(grid, scale_factor=1.02))
        self.play(FadeOut(VGroup(title, subtitle, grid, k_brace, k_lbl, s_brace, s_lbl,
                                 r_brace, r_lbl, v_note)))


class S10_Alpha01Build(VoiceoverScene):
    """The Series III generating array that produces the Part 2 design."""

    def construct(self):
        voice(self)
        title = Text("Part 2  -  the generating array", font_size=34,
                     weight="BOLD").to_edge(UP)
        subtitle = Text("Series III   (Patterson & Williams 1976)",
                        font_size=26, color=C_HILITE).next_to(title, DOWN, buff=0.15)
        self.add(title, subtitle)
        agrid = NumberGrid(D2["alpha"], color=C_SEED).scale(1.2)
        with self.voiceover(
            text="How is that design built? For this case -- three replications and an even "
            "number of blocks -- the paper gives a Series three construction, the same family "
            "FieldHub uses."
        ):
            self.play(Create(agrid))
        with self.voiceover(
            text="Here is its generating array: five rows, one per plot in a block, and three "
            "columns, one per replication. Developing it modulo six, and shifting the rows by "
            "zero, five, ten, fifteen and twenty, produces exactly the thirty-variety design "
            "we just saw."
        ):
            self.play(Indicate(agrid))
        self.play(FadeOut(VGroup(title, subtitle, agrid)))


class S11_SingleCoOccur(VoiceoverScene):
    def construct(self):
        voice(self)
        p = D2["params"]
        r, s, k = p["r"], p["s"], p["k"]
        a, b = D2["track_pair"]
        bm = _block_matrix(D2)
        grid = NumberGrid(bm, color=C_SEED, cell=P2_CELL, fs=P2_FS)
        title = Text("α(0,1): every pair co-occurs at most once",
                     font_size=30).to_edge(UP)
        # Reproduce the Part 2 design layout as numbered block-cards, coloured per replication.
        for c in range(r * s):
            grid.column(c).shift(RIGHT * _split_shift(c, r, s, P2_REPGAP, P2_COLGAP))
        rep_labels = VGroup()
        for j in range(r):
            for c in range(j * s, (j + 1) * s):
                for rr in range(k):
                    grid.square(rr, c).set_stroke(REP_COLORS[j], width=3)
            grp = VGroup(*(grid.column(c) for c in range(j * s, (j + 1) * s)))
            rep_labels.add(Text(f"Rep {j + 1}", font_size=24, color=REP_COLORS[j])
                           .next_to(grp, UP, buff=0.3))
        self.add(title, grid, rep_labels)

        def col_vals(c):
            return [bm[i][c] for i in range(k)]

        both = [c for c in range(r * s) if a in col_vals(c) and b in col_vals(c)]
        boxes = VGroup()
        between = VGroup()
        for c in both:
            ra, rb = col_vals(c).index(a), col_vals(c).index(b)
            for row in (ra, rb):
                boxes.add(SurroundingRectangle(grid.cell(row, c), color=C_HILITE,
                                               buff=0.05, stroke_width=4))
            for row in range(min(ra, rb) + 1, max(ra, rb)):   # plots sitting between them
                between.add(grid.cell(row, c))
        with self.voiceover(
            text=f"Track varieties {a} and {b}. They share exactly one block -- here -- "
            "and they never appear together anywhere else."
        ):
            self.play(Create(boxes))
        with self.voiceover(
            text="And notice they are not neighbours: another plot sits between them in the "
            "block. Co-occurring is about sharing a block, not about being next to each other."
        ):
            self.play(Indicate(between, color=C_OFFSET))

        counts = D2["cooccurrence_offdiag_counts"]
        n0 = counts.get("0", "0")
        n1 = counts.get("1", "0")
        n2 = counts.get("2", "0")
        total = sum(int(x) for x in counts.values())
        tally = VGroup(
            Text(f"{total} variety pairs", font_size=26, weight="BOLD"),
            Text(f"{n0} never share a block", font_size=24, color=BLUE),
            Text(f"{n1} share exactly one block", font_size=24, color=GREEN),
            Text(f"{n2} share two blocks", font_size=24, color=RED),
        ).arrange(DOWN, aligned_edge=LEFT, buff=0.16).to_edge(DOWN)
        with self.voiceover(
            text=f"Across all {total} pairs, {n0} never meet, {n1} meet in exactly one "
            f"block, and {n2} meet twice. No pair is forced to co-occur twice, because the "
            "block size does not exceed the number of blocks. That is an alpha, zero one, "
            "design."
        ):
            self.play(FadeIn(tally, shift=UP * 0.2))
        self.wait(0.3)


# ===========================================================================
class S12_CoocPart2(VoiceoverScene):
    """Co-occurrence heat-map + histogram for the Part 2 v=30 α(0,1) design."""

    def construct(self):
        voice(self)
        p = D2["params"]
        title = Text(f"Co-occurrence structure  ({p['v']} varieties, k={p['k']}, s={p['s']})",
                     font_size=30).to_edge(UP)
        self.add(title)

        heat = make_heatmap(D2["cooccurrence_matrix"], cell=0.18)
        heat.to_edge(LEFT, buff=0.9).shift(DOWN * 0.3)
        legend = make_legend().next_to(heat, DOWN, buff=0.3)
        with self.voiceover(
            text="Here is the co-occurrence heat map for this alpha zero one design. Thirty "
            "varieties -- blue and green only. There is no red, because no pair shares two "
            "blocks."
        ):
            self.play(FadeIn(heat), FadeIn(legend))

        counts = D2["cooccurrence_offdiag_counts"]
        total = sum(int(x) for x in counts.values())
        hist = make_histogram(counts).scale(0.9)
        hist.to_edge(RIGHT, buff=1.0).shift(DOWN * 0.3)
        acls = alpha_label("α(0, 1) design").next_to(hist, DOWN, buff=0.45)
        with self.voiceover(
            text=f"Of the {total} pairs, {counts['0']} never meet, {counts['1']} meet once, "
            "and none meet twice."
        ):
            self.play(FadeIn(hist, shift=UP * 0.2))
            self.play(Write(acls))
        self.wait(0.4)


class S13_RealTrial(VoiceoverScene):
    def construct(self):
        voice(self)
        rt = D2["real_trial"]
        p = rt["params"]
        cc = rt["cooccurrence_offdiag_counts"]
        npairs = rt["n_pairs"]
        title = Text("Your trial, as an α(0,1) design",
                     font_size=34, weight="BOLD").to_edge(UP)
        params = VGroup(
            Text(f"{rt['software']}  -  Series {rt['series']} construction", font_size=26),
            Text(f"{p['v']} entries    k = {p['k']}    s = {p['s']} blocks/rep    "
                 f"r = {p['r']} reps", font_size=26),
        ).arrange(DOWN, buff=0.25).next_to(title, DOWN, buff=0.6)
        stats = VGroup(
            Text(f"{npairs} variety pairs", font_size=28, weight="BOLD"),
            Text(f"{cc['0']} never share a block", font_size=26, color=BLUE),
            Text(f"{cc['1']} share exactly one block", font_size=26, color=GREEN),
            Text(f"{cc['2']} share two blocks", font_size=26, color=RED),
        ).arrange(DOWN, aligned_edge=LEFT, buff=0.2).next_to(params, DOWN, buff=0.55)

        with self.voiceover(
            text="This scales straight to your real trial. Fifty entries, blocks of five, "
            "ten blocks per replication, three replications, generated by FieldHub with "
            "the same Series three construction."
        ):
            self.play(Write(title))
            self.play(FadeIn(params, shift=UP * 0.2))
        with self.voiceover(
            text=f"Because the block size, five, is at most the number of blocks, ten, it "
            f"is an alpha zero one design. Of all {npairs} pairs, {cc['0']} never meet, "
            f"{cc['1']} meet exactly once, and not one pair is forced to co-occur twice."
        ):
            self.play(FadeIn(stats, shift=UP * 0.2))
        with self.voiceover(
            text="Same recipe, same randomization, larger field. Construct once, randomize "
            "per site."
        ):
            self.play(Indicate(stats[3]))
        self.wait(0.5)


# ===========================================================================
class S14_CoocReal(VoiceoverScene):
    """Co-occurrence heat-map + histogram for the real v=50 FieldHub α(0,1) design."""

    def construct(self):
        voice(self)
        rt = D2["real_trial"]
        title = Text("Co-occurrence structure  (your 50-entry FieldHub design)",
                     font_size=30).to_edge(UP)
        self.add(title)

        heat = make_heatmap(rt["cooccurrence_matrix"], cell=0.12)
        heat.to_edge(LEFT, buff=0.8).shift(DOWN * 0.3)
        legend = make_legend().next_to(heat, DOWN, buff=0.3)
        with self.voiceover(
            text="The same heat map for your FieldHub design. Fifty entries, so twenty-five "
            "hundred cells -- and only blue and green. There is no red anywhere."
        ):
            self.play(FadeIn(heat), FadeIn(legend))

        hist = make_histogram(rt["cooccurrence_offdiag_counts"]).scale(0.9)
        hist.to_edge(RIGHT, buff=1.0).shift(DOWN * 0.3)
        with self.voiceover(
            text="Nine hundred twenty-five pairs never meet, three hundred meet once, and "
            "not a single pair meets twice."
        ):
            self.play(FadeIn(hist, shift=UP * 0.2))
        acls = alpha_label("α(0, 1) design").next_to(hist, DOWN, buff=0.45)
        with self.voiceover(text="A clean alpha, zero one, design."):
            self.play(Write(acls))
        self.wait(0.4)
