# Storyboard & Narration — α-Lattice Design Randomization

Worked example throughout: **v = 20 varieties, r = 3 replications, s = 4 blocks/rep, k = 5 plots/block** (Patterson & Williams 1976, Table 1).

Narration uses the macOS `say` voice via `manim-voiceover` (local `say_service.py`). Each scene = one Manim `Scene`; the `with self.voiceover(text=...)` block wraps the animations so they sync to the spoken line. Scenes are **S01–S14**, in playback order (matching the class names in `alpha_lattice.py`). Keep on-screen numbers identical to `design_data*.json`.

---

## Scene 1 — What an α-design optimizes (introduction)  (`S01_Optimize`)
**On screen:** Title "What an alpha-design optimizes". A 10-plot field strip with an opacity gradient (a fertility trend). First a gold box around the whole strip ("one big block: large error variance"); then it splits into 5 small blue boxes ("small blocks: nearly uniform soil, small variance"). Cut to the formula `precision of a comparison ∝ (1/σ²) × E`, with two notes: "block size sets the variance (experimenter's choice)" and the highlighted "E = fraction of precision kept, 0 < E ≤ 1".

**Narration:**
> *(1)* Why split a trial into small incomplete blocks at all? It comes down to one trade-off. Picture a field with a fertility trend running across it.
> *(2)* If a whole replication is one big block, plots at opposite ends sit on very different soil, so the background variance is large, and real variety differences are hard to see.
> *(3)* Break it into small blocks and each block covers nearly uniform soil. The variance drops, and comparisons within a block get sharper.
> *(4)* But there is a catch. Varieties in different blocks are no longer compared directly, so some information is confounded with block differences. The precision of a comparison depends on two things: the variance σ², and an efficiency factor E.
> *(5)* The block size sets σ² — the experimenter's choice. The efficiency factor E, between 0 and 1, is the fraction of comparison precision you keep despite the blocking.
> *(6)* The α-design has a single objective: **maximize efficiency** — the factor E, the precision of variety comparisons kept despite the blocking. It reaches that goal in just one way: **by minimizing co-occurrences** — keeping any two varieties from sharing a block more often than they must. (On-screen callout: "one goal: maximize efficiency (E) / by minimizing co-occurrences".) Central point = maximize E; the means = minimize co-occurrences (paper §5: minimize the number of concurrences greater than one).

---

## Scene 2 — The breeder's problem  (`S02_Problem`)
**On screen:** Title "Alpha-Lattice Designs" / "Patterson & Williams (1976)". Then a **preview of the final design** — colour-by-rep, split block-cards, **no numbers**: 3 reps (blue/green/gold), each split into s=4 blocks, each block a solid column of k=5 plots (no space between plots in a block; gaps between blocks and larger gaps between reps). Brace annotations span the dimensions: **k = 5** (vertical, one block), **s = 4** (over one rep), **r = 3** (across all reps), with **v = k × s = 20** below. (Replaces the old chips.)

**Narration:**
> *(1)* Suppose we want to compare twenty varieties in a field trial, with three full replications.
> *(2)* A whole replication in one long strip is noisy, because distant plots differ in soil and moisture. So we break each replication into smaller, more uniform incomplete blocks: four blocks of five plots each.
> *(3)* Three replications in all, and twenty varieties, since twenty is five times four. The challenge is choosing which varieties share a block, so comparisons stay balanced. Patterson and Williams gave a simple recipe.

---

## Scene 3 — The generating array α  (`S03_GeneratingArray`)
**On screen:** Draw the 5×3 array α, with a static parameter chip pinned at the bottom (not narrated): "k = 5   s = 4   r = 3". (No spoken recap — Scene 2 already set up 20 varieties / 3 reps / blocks of 5.) Highlight that entries are residues mod 4 (0,1,2,3). Box the first row and first column to show the "reduced" form (all zeros there).
```
0 0 0
0 1 2
0 2 3
0 3 1
0 3 2
```
**Narration:**
> The recipe starts with a generating array: five rows for the five plots in a block, three columns for the three replications. Every entry is a number from zero to three, that is, a residue modulo four, where four is the number of blocks in a replication. This particular array is in reduced form: the whole first row and first column are zero. Those few numbers are the entire seed of the design. And they are not arbitrary: the array is **precomputed** — looked up from tables, given by a formula, or found by search — and chosen so the finished design **maximizes the efficiency factor E**. Picking the array is the only clever step; everything after it is mechanical.

---

## Scene 4 — Cyclic development → α\*  (`S04_CyclicDevelopment`)
**On screen:** Take each column of α and expand it rightward into 4 columns: add 1 mod 4 each step. Animate the "+1, wrap 3→0" for one column slowly, then fast-forward the rest. Result is the 5×12 intermediate array α\*:
```
0 1 2 3 | 0 1 2 3 | 0 1 2 3
0 1 2 3 | 1 2 3 0 | 2 3 0 1
0 1 2 3 | 2 3 0 1 | 3 0 1 2
0 1 2 3 | 3 0 1 2 | 1 2 3 0
0 1 2 3 | 3 0 1 2 | 2 3 0 1
```
Use vertical dividers to mark the three groups of 4 columns.

**Narration:**
> Now we develop each column cyclically. Take a column, then make three more copies, adding one to every number each time — and whenever a three would become four, it wraps back around to zero. Watch the first column: all zeros, then ones, then twos, then threes. Each original column blossoms into a group of four. The three columns of the seed become twelve columns in this intermediate array. The wrap-around is the key trick: it spreads the varieties evenly, so no two end up paired more often than necessary.

---

## Scene 5 — Row offsets → variety labels  (`S05_RowOffsets`)
**On screen:** Show offsets entering from the left of each row: +0, +4, +8, +12, +16. Animate each row's entries incrementing to their final variety labels 0–19. End with the labelled 5×12 grid.

**Narration:**
> So far every number is still small, between zero and three. To turn them into actual variety labels from zero to nineteen, we shift each row by a fixed amount. The first row stays as is. We add four to the second row, eight to the third, twelve to the fourth, and sixteen to the fifth. Now every one of the twenty varieties appears, and — because of how we built the rows — each variety shows up exactly once in every group of four columns.

---

## Scene 6 — Blocks & replications  (`S06_Blocks`)
**On screen:** **Opens on the exact blue variety matrix that Scene 5 (S05) ended on (seamless cut).** Then a smooth transition: (1) the 12 columns **split** into three groups at the replication boundaries (gaps open), (2) each group's cells are **recoloured** to its rep colour (Rep 1 blue, Rep 2 green, Rep 3 gold), (3) **Rep 1 / 2 / 3 labels** fade in above each group. Columns = blocks, matching Table 1 (Rep I {0,4,8,12,16}…, etc.).

Then highlight varieties **4 and 16**: they share a column (block) in Rep I *and* Rep III → "co-occur twice". **They are non-adjacent within each column** (plots 8 and 12 sit between them) — teaching that co-occurrence = sharing a block, *not* neighbouring plots. Caption: "α(0, 1, 2) design".

**Narration:**
> *(split)* These twelve columns are the blocks of the design. First, watch them separate at the replication boundaries.
> *(colour)* We colour the three groups: replication one, replication two, replication three.
> *(label)* Each group of four columns is one complete replication, containing all twenty varieties. So the design is resolvable.
> *(block split)* And each column is a block — a set of five varieties grown together. *(The four columns within each rep separate into individual block-cards, gap 0.18, matching the block-card spacing of the randomization section.)*
> *(pair)* Look at varieties four and sixteen. They share a block twice — once in replication one and once in replication three. And notice they are not neighbours: other plots sit between them in the column. Co-occurring means sharing a block, not being side by side.
> *(class)* Most pairs meet zero or one times; a few meet twice. That tight co-occurrence structure makes the design efficient. This is an α(0,1,2) design.

## Scene 7 — Co-occurrence heat-map & histogram (Part 1)  (`S07_CoocPart1`)
**On screen:** A **20×20 heat-map** of the co-occurrence matrix (blue=0, green=1, red=2, grey diagonal) on the left with a colour legend; a **labelled histogram** on the right with bars 82 (blue, 0×), 96 (green, 1×), 12 (red, 2×).

**Narration:**
> *(1)* Here is the whole co-occurrence structure as a heat map. Each cell is a pair of varieties: blue if they never share a block, green if they share one, red if they share two. The grey diagonal is a variety with itself.
> *(2)* Counting them: eighty-two pairs never meet, ninety-six meet once, and twelve meet twice. Those twelve red cells are what make this an α(0,1,2) design.

---

## Scene 8 — Randomization → field layout, with pair-tracking  (`S08_Randomization`)
**On screen:** Take the systematic plan and apply four shuffles in sequence. **A maroon box stays locked on varieties 4 & 16 throughout** (a non-adjacent pair — plots sit between them), and a bottom badge reads "varieties 4 & 16 share a block: 2×" (the "2×" pulses at each step). The four steps:
1. **Permute replications** — boxes travel with their blocks.
2. **Permute blocks within each replication** — still two boxes, two shared blocks.
3. **Permute plots within each block** — 4 & 16 slide within their blocks (the gap between them changes) but stay in the same block: proximity ≠ sharing a block.
4. **Relabel varieties** (coprime map `new = (old·3+7) mod 20`: 7→8, 8→11) — boxes do **not** move (same blocks), only the numbers change; badge updates to "8 & 11 … 2×".

Final caption: "Construct once, randomize per site." This scene is the payoff for the co-occurrence idea: it shows *why* randomization is safe.

**Narration (one block per step, matching the code):**
> *(intro)* The plan is systematic, so before planting we randomize it. But watch varieties four and sixteen, which share a block twice — and note they are never side by side; other plots sit between them. The maroon boxes will follow them through every shuffle.
> *(1)* First we shuffle the order of the three replications. The boxes move with their blocks, but it is still the same two blocks.
> *(2)* Then, within each replication, we shuffle the order of the blocks. Two boxes still, two shared blocks still.
> *(3)* Within each block, we shuffle the plot positions. Seven and eight slide around inside their blocks, but they stay together.
> *(4)* Finally we relabel every variety. Four becomes nineteen, and sixteen becomes fifteen. The boxes do not move, because it is the very same two blocks. Only the names changed.
> *(close)* The co-occurrence count never changed: still exactly two. Randomization preserves the design's co-occurrence structure, while spreading any field trend randomly across varieties. A fair, efficient, ready-to-plant alpha-lattice layout.

---

# PART 2 — the α(0,1) design (k ≤ s ⇒ single co-occurrences)

Second worked example: a **different design** — **v = 30, r = 3, s = 6, k = 5** — built with the paper's **Series III** construction (r=3, even s), the *same family as the FieldHub trial* (k=5, s=10). It is **not** Part 1's trial with extra blocks (that would just need more varieties); it is a separate design used to illustrate the **improved α(0,1) lattice**. When **k ≤ s** (here 5 ≤ 6), no pair is forced to share a block twice → α(0,1). The improvement over Part 1's α(0,1,2) is *abstract*: **fewer co-occurrences**, spreading varieties more evenly (higher E). Verified data in `design_data_alpha01.json` (co-occurrence counts 255 / 180 / 0).

## Scene 9 — The α(0,1) design (Scene-2 style)  (`S09_Alpha01Design`)
**Point of Part 2:** a separate 30-variety design where **block size k ≤ number of blocks per replication s**, so every pair co-occurs at most once → the **α(0,1)** lattice. The improvement is the reduced co-occurrence (better spread / higher E), not "more blocks of the same size."

**On screen:** Title "Part 2 — the α(0,1) design", subtitle "when k ≤ s, no pair co-occurs more than once   (here k = 5, s = 6)". Show the v=30 design **as in Scene 2** — colour-by-rep, split block-cards, **no internal numbers**: 3 reps (blue/green/gold), each s=6 blocks, each block a solid column of k=5 plots. Brace annotations: **k = 5** (one block), **s = 6** (one rep), **r = 3** (all reps), with **v = k × s = 30** below. (Same block-card layout that Scene 11 then uses, with numbers.)

**Narration:**
> *(1)* Part two is a different design — a separate trial, with thirty varieties — illustrating the improved lattice you get when the block size k is at most the number of blocks s. Here, blocks of five plots and six blocks per replication: k = 5, s = 6.
> *(2)* Three replications — thirty varieties in all, since thirty is five times six. Because k ≤ s, no pair is ever forced to share a block more than once: the improvement of an α(0,1) design over an α(0,1,2).

## Scene 10 — The generating array  (`S10_Alpha01Build`)
**On screen:** Title "Part 2 — the generating array". Show the 5×3 Series III generating array `[[0,0,0],[0,1,3],[0,2,1],[0,3,4],[0,4,2]]` (with numbers; a Series III generator — the paper gives the *method*, not this exact array).

**Narration:**
> *(1)* How is that design built? For this case — three replications and an even number of blocks — the paper gives a Series III construction, the same family FieldHub uses.
> *(2)* Here is its generating array: five rows, one per plot in a block, and three columns, one per replication. Developing it modulo six, and shifting the rows by 0, 5, 10, 15 and 20, produces exactly the thirty-variety design we just saw.

## Scene 11 — Single co-occurrences  (`S11_SingleCoOccur`)
**On screen:** Repeats the Part 2 design layout (Scene 9), now as numbered block-cards grouped by rep with Rep 1/2/3 labels. Box the tracked pair (`0` & `12`) in the *one* block they share — non-adjacent (variety `6` between them, pulsed gold). Bottom tally (palette blue/green/red): "435 variety pairs / 255 never / 180 once / 0 twice".

**Narration:**
> *(1)* Track varieties zero and twelve. They share exactly one block — here — and they never appear together anywhere else.
> *(1b)* And notice they are not neighbours: another plot sits between them in the block. Co-occurring is about sharing a block, not about being next to each other.
> *(2)* Across all 435 pairs, 255 never meet, 180 meet in exactly one block, and 0 meet twice. No pair is forced to co-occur twice, because the block size does not exceed the number of blocks. That is an α(0,1) design.

## Scene 12 — Co-occurrence heat-map & histogram (Part 2)  (`S12_CoocPart2`)
**On screen:** A **30×30 heat-map** of the Part 2 co-occurrence matrix (blue=0, green=1, grey diagonal; **no red**) with legend, and a **histogram** (255 blue 0×, 180 green 1×, 0 red 2×) with an "α(0,1) design" label beneath.

**Narration:**
> *(1)* Here is the co-occurrence heat map for this α(0,1) design. Thirty varieties — blue and green only. There is no red, because no pair shares two blocks.
> *(2)* Of the 435 pairs, 255 never meet, 180 meet once, and none meet twice.

## Scene 13 — Your real trial  (`S13_RealTrial`)
**On screen:** Title "Your trial, as an α(0,1) design". Panel: "FieldHub — Series III construction; 50 entries, k=5, s=10 blocks/rep, r=3 reps". Stats: "1225 variety pairs / 925 never share a block / 300 share exactly one block / 0 share two blocks".

**Narration:**
> *(1)* This scales straight to your real trial. Fifty entries, blocks of five, ten blocks per replication, three replications, generated by FieldHub with the same Series three construction.
> *(2)* Because the block size, five, is at most the number of blocks, ten, it is an alpha zero one design. Of all twelve hundred twenty-five pairs, nine hundred twenty-five never meet, three hundred meet exactly once, and not one pair is forced to co-occur twice.
> *(3)* Same recipe, same randomization, larger field. Construct once, randomize per site.

## Scene 14 — Co-occurrence heat-map & histogram (FieldHub)  (`S14_CoocReal`)
**On screen:** A **50×50 heat-map** of the FieldHub co-occurrence matrix (blue=0, green=1, grey diagonal; **no red**) with legend, and a **histogram** with bars 925 (blue, 0×), 300 (green, 1×), 0 (red, 2×).

**Narration:**
> *(1)* The same heat map for your FieldHub design. Fifty entries, so twenty-five hundred cells — and only blue and green. There is no red anywhere.
> *(2)* Nine hundred twenty-five pairs never meet, three hundred meet once, and not a single pair meets twice. A clean α(0,1) design.

---

## Render notes
- Voiceover: macOS `say` "Evan (Enhanced)" via the custom `say_service.py` (rate 155; sentence/comma/dash/trailing pauses). Total ~8.5 min.
- `build.sh` concatenates per-scene MP4s with an audio re-encode + `aresample=async=1` (gapless, no drift).
- Heat-map palette is fixed in `alpha_lattice.py` (`HEAT`): 0=blue, 1=green, 2=red, diagonal=grey.
- All arrays/blocks/matrices are read from the `design_data*.json` files — never hard-code numbers in the scene file.
