# Maize Combining Ability & Phosphorus-Use Efficiency (Line × Tester)

📊 **Rendered analysis site (GitHub Pages): <https://faustovrz.github.io/gcamoz/>**

A structured R/Quarto repository for the analysis of a low-phosphorus maize line × tester
trial, and for teaching General/Specific Combining Ability (GCA/SCA). It supports Stélio
Boaventura Nuvunga's PhD research, **"Genetic Analysis of Tolerance to Low Phosphorus in
Maize,"** and maize breeding for low-phosphorus environments in Mozambique.

***

## Paper analysis series (PUE hybrid study)

The manuscript's analytical framework (Joe Gage's plan). Each analysis is a standalone
notebook that feeds the paper's **Methods** and its **figures/tables**; all share one data
pipeline (`R/prep_traits.R`) and the same 11 phenotypes. Start from the
[landing page](https://faustovrz.github.io/gcamoz/).

| # | Notebook | Purpose |
|---|---|---|
| 1 | `descriptive_stats.qmd` | Data quality (missing, IQR outliers), distributions, normality (Shapiro–Wilk, Q–Q) |
| 2 | `phosphorus_effects.qmd` | Optimal-vs-stress phosphorus contrast at Chokwe (the paired site) |
| 3 | `combining_ability_anova.qmd` | Line × tester ANOVA and GCA/SCA effects (superior parents & hybrids) |
| 4 | `variance_components.qmd` | REML variance components, heritability, gene action (AV, DV, H², h², BR), BLUP selection |
| 5 | `correlations_pue.qmd` | Genetic correlations and regression — traits that drive PUE (indirect selection) |
| 6 | `stability_ammi_gge.qmd` | AMMI & GGE biplots, stability (WAAS/WAASY) for GY and PUE |

**The 11 analysed phenotypes**, grouped as in the plan:

* **A – Selection:** `GY` (grain yield, t ha⁻¹), `PUE` (grain yield / applied P)
* **B – Physiological:** `HI` (harvest index), `RSR` (root:shoot), `RW` (root weight),
  `SW` (shoot dry weight), `TDM` (total dry mass = SW + RW)
* **C – Auxiliary:** `PH` (plant height), `DTA` (days to anthesis), `DTS` (days to silking),
  `ASI` (anthesis–silking interval)

Column names and units are documented in `data/multilocation_dictionary.csv`; derived traits
(DTA, DTS, ASI) and design factors are built once in `R/prep_traits.R`.

***

## Breeding and experimental design

### Mating scheme — 7 × 7 Line × Tester (North Carolina Design II)

* **Females (lines):** 7 CIMMYT low-phosphorus donor lines — `CML364`, `CML366`, `CML434`,
  `CML435`, `CML439`, `CML530`, `CML532`.
* **Males (testers):** 7 IIAM donor testers, named by their **IIAM IDs** — `EN17`, `EN20`,
  `EN21`, `EN25`, `EN31`, `EN32`, `EN64` (legacy NCDII IDs `MOZL3`/`MOZL5`/`MOZL6`/`MOZL7`/
  `MOZL8`/`MOZL9`/`MOZL10`; full pedigrees in `multilocation_alpha_lattice`).
* **Progeny:** 49 F₁ hybrids (no selfs/reciprocals) + a replicated check, `NAMULONGUE`
  (`gen = 0`).

### Field trial & phosphorus treatments

`data/multilocation.csv` is a balanced trial analysed as an **RCBD with row/column spatial
control** (no repeated incomplete blocks were detectable — see `block_assignment_diagnostics`).

* **5 environments** across **4 physical locations** (`loc`): Chokwe was grown at **two
  phosphorus levels**, giving `env` = location × treatment:
  * **Optimal (50 kg P ha⁻¹):** `CHOKWE OPT`, `MANIQUENIQ OPT`
  * **Stress (10 kg P ha⁻¹):** `CHOKWE STS`, `NHACOONGO STS`, `SUSSUNDENGA STS`
  * The applied P is recoverable as `GY × 1000 / PUE` (exactly 10 or 50); the P effect is
    estimated only at Chokwe, the one paired site.
* **3 replications** per environment (`I`, `II`, `III`); **50 entries** each (49 hybrids +
  check), sown in a **5 × 10 serpentine** layout (field row/col recovered from planting order).
* **750 plots total** (49 hybrids × 15 + 15 checks) — balanced, complete 7 × 7 factorial.

***

## Earlier / teaching notebooks

* `gca_calculation.qmd` — teaching GCA/SCA on the **simulated** dataset
  `data/maize_factorial_yield.csv` (Line × Tester under RCBD, 2 environments).
* `multilocation_alpha_lattice.qmd` — the real trial's RCBD field layout + full line × tester
  GCA/SCA, variance components and heritability on grain yield.
* `rcbd_vs_alpha_lattice.qmd`, `block_assignment_diagnostics.qmd` — design diagnostics
  (RCBD vs alpha-lattice; co-occurrence and permutation tests).
* `male_gca_trait_sweep.qmd`, `male_gca_singularity_and_variance_routing.md` — trait-by-trait
  tester-GCA refit and the germplasm interpretation (why male GCA is ~null for yield).
* `pue_distribution.qmd` — the PUE bimodality check (two P treatments).
* `diallel_gca_sca_analysis.qmd`, `iris_anova_cld.qmd`, `retrieve_CML_genesys.qmd` — reference
  notebooks (full diallel, ANOVA + compact-letter display, CIMMYT CML passport retrieval).

***

## Repository structure

```text
gcamoz/
├── README.md                       # this file
├── _quarto.yml                     # Quarto project (output-dir: docs)
├── index.qmd                       # landing page
├── R/
│   └── prep_traits.R               # shared data prep: factors, derived traits, 11-trait dictionary
├── data/
│   ├── multilocation.csv           # cleaned real data (750 rows, abbreviated schema)
│   ├── multilocation_dictionary.csv# column dictionary (name, unit, role)
│   ├── MULTILOCATION DATAS GENERAL STANDARD.xlsx  # source workbook from Stelio
│   └── maize_factorial_yield.csv   # simulated teaching dataset
├── descriptive_stats.qmd           # paper NB1
├── phosphorus_effects.qmd          # paper NB2
├── combining_ability_anova.qmd     # paper NB3
├── variance_components.qmd         # paper NB4
├── correlations_pue.qmd            # paper NB5
├── stability_ammi_gge.qmd          # paper NB6
├── multilocation_alpha_lattice.qmd # + rcbd_vs_alpha_lattice, block_assignment_diagnostics,
│                                    #   male_gca_trait_sweep, pue_distribution, gca_calculation, …
└── docs/                           # rendered HTML (served by GitHub Pages)
```

## Required R packages

```r
install.packages(c("tidyverse", "sommer", "lme4", "lmerTest", "emmeans",
                   "Hmisc", "metan", "ggbeeswarm", "knitr", "kableExtra"))
```

* `sommer`, `lme4`/`lmerTest`, `emmeans` — mixed models (variance components, BLUPs, treatment tests)
* `Hmisc` — correlation matrices with p-values; `metan` — AMMI/GGE stability
* `ggbeeswarm` — beeswarm + transparent-boxplot trait figures
* `tidyverse`, `knitr`, `kableExtra` — data wrangling, tables

## Compilation

This is a Quarto project (`_quarto.yml`), with `docs/` as the output directory that GitHub
Pages serves. Render everything, or a single notebook:

```bash
quarto render                       # all registered notebooks -> docs/
quarto render variance_components.qmd
```

***

## References

Bibliographic metadata lives in [`references.bib`](references.bib).

* Comstock, R. E., & Robinson, H. F. (1948). The Components of Genetic Variance in
  Populations of Biparental Progenies. *Biometrics*, 4(4), 254–266.
  [doi:10.2307/3001412](https://doi.org/10.2307/3001412)
* Isik, F., Holland, J., & Maltecca, C. (2017). *Genetic Data Analysis for Plant and Animal
  Breeding*. Springer. [doi:10.1007/978-3-319-55177-7](https://doi.org/10.1007/978-3-319-55177-7)
