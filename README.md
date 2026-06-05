# Maize Combining Ability Field-Trial Analysis (Line × Tester GCA/SCA)

A structured R/Quarto repository designed to instruct plant-breeding students in maize field-trial analysis, General Combining Ability (GCA), and Specific Combining Ability (SCA).

This project is tailored specifically to Stélio Boaventura Nuvunga's PhD research on **"Genetic Analysis of Tolerance to Low Phosphorus in Maize"**.

***

## Two connected workflows

1.  **Simulated teaching workflow** — `gca_calculation.qmd` on `data/maize_factorial_yield.csv`: a 7 × 7 Line × Tester / North Carolina Design II dataset under an **RCBD** assumption, used to teach the GCA/SCA calculations step by step. It uses placeholder parent IDs (`MOZ1`–`MOZ7`, and a `CML…440…` female set) and **2 environments** (`env1` normal-P, `env2` low-P stress), 294 plot rows.

2.  **Stelio's real multilocation workflow** — `multilocation_alpha_lattice.qmd` on `data/multilocation.csv`: the real plot-level field trial, which now carries the line × tester **pedigree** and a replicated **check**, analyzed as a **randomized complete block design (RCBD) with row/column spatial control**.

***

## Breeding and Experimental Design (real trial)

### Mating scheme

A **7 × 7 Factorial Cross (Line × Tester / North Carolina Design II)**:

*   **Females (lines):** 7 CIMMYT low-phosphorus donor lines — `CML364`, `CML366`, `CML434`, `CML435`, `CML439`, `CML530`, `CML532`.

*   **Males (testers):** 7 Mozambican local testers — `MOZL3`, `MOZL5`, `MOZL6`, `MOZL7`, `MOZL8`, `MOZL9`, `MOZL10`.

*   **Progeny:** 49 $F_1$ hybrid families (no selfs or reciprocals), plus a replicated commercial check, **`NAMULONGUE`** (`GEN = 0`).

The real file now carries the pedigree directly in the `Female`, `Male`, and `Female x Male` columns, so **no external genotype-to-cross key is needed**. (Note: the *simulated teaching* dataset uses the placeholder IDs `CML…440…` and `MOZ1`–`MOZ7`; an earlier draft mistakenly listed `CML440` for the real set — the real female set includes `CML530`.)

### Field-trial layout

`data/multilocation.csv` is a balanced multilocation field trial analyzed as an **RCBD**:

*   **50 entries per replication:** 49 hybrids (`GEN = 1`–`49`) + 1 `NAMULONGUE` check (`GEN = 0`).

*   **4 environments/locations:** `CHOKWE STS`, `CKOKWE OPT`, `NHACOONGO`, `SUSSUNDENGA`.

*   **3 replications** per environment: `I`, `II`, `III`.

*   **5 × 10 serpentine layout:** each replication is sown as **5 columns × 10 rows** (50 plots), first plot top-left, snaking down along the `ENTRY` planting order. Each plot is a **full field row of maize** (not a pot). Field **row** (1–10) and **column** (1–5) are recovered from the planting order and fit as nested spatial terms.

*   **600 plots total** (49 hybrids × 12 + 12 checks) — balanced, no duplicates. *(An earlier file version had 597 rows with the check mislabeled as a duplicate `GEN = 49`; the cleaned file fixes this with the check as its own `GEN = 0` entry.)*

**Design status — alpha-lattice vs. RCBD.** Stelio describes the trial as an **alpha (incomplete-block) design**, but the data show **no incomplete, repeated sub-blocks of hybrids** — the hybrids appear fully randomized within each 5 × 10 replication — so it is analyzed here as an **RCBD with row/column spatial control**. If Stelio provides the incomplete-block (sub-block) assignment per plot, an `env:rep:block` term can be added and compared (logLik / AIC) against the row/column model.

### Evaluated traits

Grain yield (`Yield t/ha`) plus `yield g/plant (g)`, `Plant Number`, `Ear weight (g)`, `grain weight (g)`, `Plant Height (m)`, `Culm Diameter (cm)`, flowering dates (`F. D. Female`, `F. D. Male`), `Dry matter mass (g)`, and `Root Weight (g)`.

***

## Repository Structure

```text
gcamoz/
├── README.md                          # Repository documentation (this file)
├── _quarto.yml                        # Quarto project config (output-dir: docs)
├── references.bib                     # Bibliography (Comstock 1948; Isik et al. 2017)
├── iris_anova_cld.qmd                 # Reference: ANOVA and Compact Letter Display
├── diallel_gca_sca_analysis.qmd       # Reference: 7x7 full-diallel combining ability
├── gca_simulation.qmd                 # Build and export the simulated teaching dataset
├── gca_calculation.qmd                # Teaching GCA/SCA on the simulated dataset (RCBD, 2 envs)
├── multilocation_alpha_lattice.qmd    # Stelio's real RCBD multilocation line x tester analysis
├── retrieve_CML_genesys.qmd           # Retrieve CIMMYT CML passport data / IDs from Genesys
├── R/
│   ├── genesys_auth.R                 # Google-login token bridge for the Genesys API
│   └── retrieve_cml.R                 # CML retrieval helpers
├── data/
│   ├── maize_factorial_yield.csv      # Simulated dataset (294 plot rows)
│   └── multilocation.csv              # Stelio's real plot-level data (600 rows)
├── output/
│   ├── cml_ids.csv                    # Retrieved CIMMYT CML IDs
│   └── genesys_passport_raw.rds       # Raw Genesys passport data
└── docs/                              # Rendered HTML (also served by GitHub Pages)
    ├── iris_anova_cld.html
    ├── diallel_gca_sca_analysis.html
    ├── gca_simulation.html
    ├── gca_calculation.html
    ├── multilocation_alpha_lattice.html
    └── retrieve_CML_genesys.html
```

### Key notebooks

*   `gca_calculation.qmd` & `docs/gca_calculation.html`: Teaching notebook on the **simulated** dataset (Line × Tester under an RCBD assumption, 2 environments). Walks through Exploratory Data Analysis, fixed-effects arithmetic GCA/SCA, a `sommer` mixed model, Breeding Values ($BV$) and Genotypic Values ($GV$), and heritability (narrow-sense, broad-sense, family-mean, and an across-environment G × E model).

*   `multilocation_alpha_lattice.qmd` & `docs/multilocation_alpha_lattice.html`: Stelio's **real** data. Documents the RCBD **5 × 10 serpentine** field layout, then runs the full line × tester pipeline on `Yield t/ha` across all four environments:

    *   *EDA:* the 7 × 7 mean-yield grid and per-environment yield distributions (with the check overlaid).

    *   *Fixed-effects arithmetic GCA/SCA* from family means.

    *   *Mixed-effects model* (`sommer::mmer`): `env` fixed; `env:rep`, `env:rep:row`, `env:rep:col` design terms; `female`, `male`, `cross` and their `:env` interactions — yielding GCA/SCA BLUPs and variance components.

    *   *Breeding and Genotypic Values* for ranking hybrids.

    *   *Entry-mean heritability with G × E*, plus a step-by-step checklist for Stelio.

*   `retrieve_CML_genesys.qmd` + `R/`: retrieve CIMMYT CML founder passport data / accession IDs from the [Genesys](https://www.genesys-pgr.org/) plant genetic resources database (with a Google-login token bridge in `R/genesys_auth.R`).

*   `gca_simulation.qmd`, `diallel_gca_sca_analysis.qmd`, `iris_anova_cld.qmd`: supporting/reference notebooks (data simulation, full-diallel combining ability, and an ANOVA + compact-letter-display example).

### Data files

*   `data/maize_factorial_yield.csv`: simulated grain yields, 294 rows, columns `env`, `block`, `female`, `male`, `cross`, `yield`. Encodes inbreeding depression, deterministic GCA effects, specific heterosis (SCA), and GCA × Environment interactions (amplified CIMMYT tolerance under low-P stress).

*   `data/multilocation.csv`: Stelio's real plot-level file, **600 rows**, columns `ENTRY`, `ENV`, `REP`, `GEN`, `Female`, `Male`, `Female x Male`, `Yield t/ha`, `yield g/plant (g)`, `Plant Number`, `Ear weight (g)`, `grain weight (g)`, `Plant Height (m)`, `Culm Diameter (cm)`, `F. D. Female`, `F. D. Male`, `Dry matter mass (g)`, `Root Weight (g)`. `ENTRY` is the field planting order.

***

## Required R Packages

```R
install.packages(c("tidyverse", "sommer", "knitr", "kableExtra"))
```

*   `tidyverse` — data manipulation, reshaping, and plotting.

*   `sommer` — quantitative-genetics mixed-model package; the `mmer` solver partitions variance components and computes BLUPs (used in `gca_calculation.qmd`, `multilocation_alpha_lattice.qmd`, and `diallel_gca_sca_analysis.qmd`).

*   `knitr` and `kableExtra` — render styled tables.

*(The Genesys retrieval notebook `retrieve_CML_genesys.qmd` has its own HTTP/JSON dependencies — see that notebook and `R/`.)*

***

## Compilation Instructions

This repository is a Quarto project (see `_quarto.yml`), with `docs/` as the output directory and the notebooks registered for rendering. The `docs/` folder is what GitHub Pages serves, so a fresh render publishes the site. To compile every notebook into `docs/`:

```bash
quarto render
```

To render a single notebook (the project-level `output-dir: docs` still applies):

```bash
quarto render multilocation_alpha_lattice.qmd
```

*(Alternatively, open the `.qmd` files in RStudio or Positron and run the chunks interactively.)*

***

## References

Bibliographic metadata lives in [`references.bib`](references.bib) and is rendered automatically by Quarto in both `gca_calculation.qmd` and `multilocation_alpha_lattice.qmd` (see their Phenotypic Model sections).

*   Comstock, R. E., & Robinson, H. F. (1948). The Components of Genetic Variance in Populations of Biparental Progenies and Their Use in Estimating the Average Degree of Dominance. *Biometrics*, 4(4), 254–266. [doi:10.2307/3001412](https://doi.org/10.2307/3001412)

*   Isik, F., Holland, J., & Maltecca, C. (2017). *Genetic Data Analysis for Plant and Animal Breeding*. Springer International Publishing. [doi:10.1007/978-3-319-55177-7](https://doi.org/10.1007/978-3-319-55177-7)
