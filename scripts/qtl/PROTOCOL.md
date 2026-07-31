# F₂ genetic map — streamlined protocol

**Population.** Single biparental F₂, DArTseq order `DMz26-3123`, 186 genotyped plants
(166 also phenotyped). Female parent **CML530**, male one of the seven `EN` testers —
established by the 7 × 7 line × tester design (`README.md:44-48`), not assumed. The
specific tester is still unknown and is needed only for the sign of allele effects.

**Coordinate system.** **CML530** (`Zm-CML530-REFERENCE-HiLo-1.0`), not B73. B73 is not a
parent of this cross, so its coordinates can place a marker on the wrong chromosome —
they disagree for 37 markers. `est.map` never reorders; it only estimates distances along
the order it is given, so a wrong order cannot be recovered downstream.

---

## Pipeline

Run in order. Each script writes a matching `.log` in `agent/`.

| # | script | does | key output |
|---|---|---|---|
| 1 | `phenotype_dictionary_unify.R` | trait dictionary | `data/phenotype_dictionary.csv` |
| 2 | `verify_dart_codes.R` | proves DArT code 2 = heterozygote | — |
| 3 | `data_audit.R` | inventory, marker + sample QC | — |
| 4 | `structure_check.R` | confirms a single biparental F₂ | — |
| 5 | `liftover_v4_to_v5.R` | AGPv4 → v5; writes the tag FASTA | `markers_v5.tsv`, `tags_for_cml530.fa` |
| 6 | `install_minibwa_from_source.sh` | Homebrew's build omits the GPL `-l` | — |
| 7 | `index_genomes.sh` | ~45 min, ~9.5 GB, peak 7.2 GB RSS | minibwa indexes |
| 8 | `cml530_alleles_direct.R` | CML530 allele per marker from its assembly | `CML530_marker_alleles_direct.tsv`, `tags_vs_CML530.bam` |
| 9 | `encode_abh.R` | filters + A/H/B encoding, **CML530 order** | `rqtl_gen_abh_all.csv`, `abh_all_marker_info.tsv` |
| 10 | `encode_phenotypes.R` | phenotypes in csvsr | `rqtl_phe.csv` |
| 11 | `build_map.R` | the map | `rqtl_cross_map_teonamqc.rds`, `data/f2_genetic_map.tsv` |

---

## Parameters, and why

### `encode_abh.R`

| param | value | basis |
|---|---|---|
| `MAPQ_MIN` | **30** + no `XA` | Tags that align equally well to several places get their CML530 allele read off an arbitrary copy. Markers anti-phase to both neighbours multi-map at 70.4% vs a 6.2% baseline (median MAPQ 4 vs 60, Fisher p = 2.3e-8). Applying the filter cut those from 27 to 5. |
| `MAF_MIN` | **0.15** | Hard cut is correct for an F₂: expected frequencies 0.25/0.50/0.25 are symmetric, MAF ≈ 0.5 at every informative marker. (A relative rule is needed only for skewed designs like BC1S4.) |
| `FDR_LEVEL` | **0.05** (BH) | Replaces a raw p > 0.01 cut, which had no stated error rate. BH lands at raw p = 0.0087 and removes 369 markers vs 376 before — the old value was already about right, this states what it controls. |
| `LINK_MIN` | **0.30** over K = 5 | Markers correlating with none of their neighbours contradict the encoding. Count is insensitive from 0.25–0.50; unreachable by any DArT statistic (RepAvg is flat, CallRate misses 80 of 117). |

Marker order is `(cml_chr, ref_pos)` from `tags_vs_CML530.bam`.

### `build_map.R`

| param | value | basis |
|---|---|---|
| drop inverted | markers anti-correlated with **both** neighbours | An inverted marker makes every individual read as a double recombinant, so rf → 0.5 both sides. Dropped rather than flipped — flipping assumes only the label is wrong. |
| distortion | per-chromosome `renorm_z`, z > 1.96 | FVRZ's relative rule, verbatim from TeoNAM. |
| `PRUNE_R` | **0.95** (= 2.565 cM) | LD prune replacing `findDupMarkers`, which needed identical genotypes and found only 28 markers. Teng & Xu 2026 Eq. 9: in an Fₜ, `r = 1 − 2θ`, so with Haldane `d = −50·ln(r)` — a genetic axis from correlation alone, which is what makes pruning *before* `est.map` possible. Thinning is Jena et al. 2018 MDdIS: on a 1-D axis it collapses to the exact O(n) interval-scheduling greedy. |
| `ERROR_PROB` | **0.01** | Both `est.map` rounds **and** `calc.genoprob` at scan time. The two criteria disagree — median-gap calibration points at 5e-3, the 1348–1596 target band at 2e-2 — and 1% is what the map-inflation arithmetic implies independently. Immaterial at scan time (r = 0.9978 vs 2.5e-3). |
| `MAP_FUN` | haldane | |
| `FINE_THR` | 10 | **Inert** — the singleton pass flags nothing at ≥ 10. Only 2 or 5 bite. |
| `ISLAND_GAP_CM` | **20** | 10 left a 55.4 cM hole on chr5 where `find_quirky` had removed three markers that were filling it. 20 is the only setting that clears the > 25 cM gaps without collapsing coverage. |
| `ISLAND_MAX_N` | 5 | |

`PROTECT_MARKERS="m1,m2"` in the environment exempts named markers from the distortion
filter, the prune and `find_quirky` — a diagnostic switch only; a map built with it set is
not a QC'd map.

---

## Attrition

```
5,208  SNPs in the DArT report
3,884  v4-anchored to chr 1-10
3,828  lifted 1:1 to v5, same chromosome
3,728  CML530 allele defined
2,994  MAPQ >= 30, no XA, mapped
2,991  on CML530 chr 1-10
2,067  MAF > 0.15
1,698  1:2:1 at BH FDR 0.05
1,692  one marker per 69 bp tag (CloneID)
1,688  |cor| >= 0.30 with a neighbour      -> encoded
1,683  not anti-phase to both neighbours
1,633  passes per-chromosome distortion
1,169  after the r = 0.95 LD prune
1,152  after find_quirky                   -> THE MAP
```

Drop rate 31.8% from the encoded set. Above the 10% guardrail; the prune accounts for
most of it and is deliberate.

## The map

**1,152 markers, 1,858.4 cM, 96.7% physical coverage** (2,109 of 2,180 Mb).
Median gap 0.834 cM, max 29.0 cM, one interval over 25 cM.
Marey maps monotonic and sigmoid on all ten chromosomes.

Published as `data/f2_genetic_map.tsv`.

## Known limitations

- **Length.** 1,858 cM against a 1,348–1,596 expectation, ~20% long. Affects cM positions
  and cM-width support intervals; **not** physical (Mb) intervals, which is what a
  candidate-gene search uses.
- **chr6** covers 78.4%, losing 34.3 Mb off the start — a QTL there cannot be found.
- **Power, not the map, is binding.** n = 166 phenotyped gives ~80% power only for QTL
  above 8–10% of variance. The two most heritable diallel traits (DTA/DTS, h² ≈ 0.60)
  were never scored on these plants.
