# Maize Combining Ability and Alpha-Lattice Trial Analysis

A structured R/Quarto repository designed to instruct plant breeding students in maize field-trial analysis, General Combining Ability (GCA), and Specific Combining Ability (SCA).

This project is tailored specifically to Stélio Boaventura Nuvunga's PhD research on **"Genetic Analysis of Tolerance to Low Phosphorus in Maize"**.

***

## Breeding and Experimental Design

The repository now contains two connected workflows:

1.  **Simulated teaching workflow:** a 7 × 7 Line × Tester / North Carolina Design II dataset used to explain GCA and SCA calculations step by step.

2.  **Stelio's real multilocation workflow:** `data/multilocation.csv`, a plot-level file sent by Stelio that is better represented as an **alpha-lattice field design**, not an RCBD.

### Breeding Mating Scheme

The intended genetic material remains a **7 × 7 Factorial Cross (Line × Tester / North Carolina Design II)**.

*   **Females (Lines):** 7 Low-Phosphorus tolerant inbred donor lines sourced from CIMMYT (`CML364`, `CML366`, `CML434`, `CML435`, `CML439`, `CML440`, `CML532`).

*   **Males (Testers):** 7 elite local Mozambican inbred lines (`MOZ1` to `MOZ7`) adapted to local agroecological zones.

*   **Evaluated Progeny:** $F_1$ factorial hybrids ($N = 49$ unique families). No selfs or reciprocals are represented in the primary line-by-tester analysis.

The real `multilocation.csv` file currently identifies entries only with `GEN` codes (`1` to `49`). A genotype-to-cross key is still needed before estimating GCA and SCA from the real data.

### Field Trial Layout

The simulated dataset in `data/maize_factorial_yield.csv` uses an RCBD assumption for teaching arithmetic GCA/SCA.

Stelio's real dataset in `data/multilocation.csv` fits an **alpha-lattice** interpretation:

*   **49 genotypes** (`GEN = 1` to `49`).

*   **4 environments/locations:** `CHOKWE STS`, `CKOKWE OPT`, `NHACOONGO`, and `SUSSUNDENGA`.

*   **3 replications** per environment: `I`, `II`, and `III`.

*   **Resolvable alpha lattice:** each replication is expected to contain **7 incomplete blocks × 7 plots = 49 plots**.

*   **Expected clean plot count:** $4 \times 3 \times 7 \times 7 = 588$ plots.

The CSV has **597 rows**, because `GEN = 49` is duplicated once in each replication for `CHOKWE STS`, `CKOKWE OPT`, and `SUSSUNDENGA`. `NHACOONGO` has the clean expected structure of 49 plots per replication. The duplicate records should be checked against the original fieldbook before final inference.

### Evaluated Traits

The real multilocation file includes grain yield (`Yield t/ha`) plus yield per plant, plant number, ear weight, grain weight, plant height, culm diameter, flowering dates, dry matter mass, and root weight.

***

## Repository Structure

The directory is clean, neat, and organized into data, rendered docs, and source notebooks:

```text
gca/
├── README.md                          # Repository documentation (this file)
├── iris_anova_cld.qmd                 # Reference example: ANOVA and Compact Letter Display
├── diallel_gca_sca_analysis.qmd       # Reference example: 7x7 full-diallel combining ability
├── gca_simulation.qmd                 # Source notebook: Setup of simulation parameters and export
├── gca_calculation.qmd                # Source notebook: Theoretical design, math equations, and R GCA/SCA analysis
├── multilocation_alpha_lattice.qmd    # Source notebook: Stelio's real multilocation alpha-lattice design
├── data/
│   ├── maize_factorial_yield.csv      # Exported simulated dataset (294 plot observations)
│   └── multilocation.csv              # Stelio's real multilocation plot-level data
└── docs/
    ├── iris_anova_cld.html            # Compiled HTML report for Iris ANOVA reference
    ├── diallel_gca_sca_analysis.html  # Compiled HTML report for Diallel GCA/SCA reference
    ├── gca_simulation.html            # Compiled HTML report for Maize Factorial simulation
    ├── gca_calculation.html           # Compiled HTML report for Maize Factorial GCA/SCA calculations
    └── multilocation_alpha_lattice.html # Compiled HTML report for real alpha-lattice data
```

### File Details

*   `data/maize_factorial_yield.csv`: Simulated grain yields containing 294 rows with columns: `env`, `block`, `female`, `male`, `cross`, and `yield`. The dataset incorporates inbreeding depression on the parental lines, deterministic GCA effects, specific heterosis (SCA), and GCA × Environment interactions (amplified CIMMYT tolerance under low phosphorus stress).

*   `data/multilocation.csv`: Stelio's real plot-level multilocation file. It contains 49 genotype codes across four environments and three replications, with traits including `Yield t/ha`, `yield g/plant (g)`, `Plant Number`, `Ear weight (g)`, `grain weight (g)`, `Plant Higher (m)`, `Colm Diameter (cm)`, flowering dates, dry matter mass, and root weight.

*   `gca_simulation.qmd` & `docs/gca_simulation.html`: The background script that sets up the parameters, builds the linear simulation model, and exports the data. Kept separate to keep Stelio's focus on calculations.

*   `gca_calculation.qmd` & `docs/gca_calculation.html`: Teaching notebook based on the simulated dataset. It outlines the **Breeding Scheme**, **Experimental Design**, and **Phenotype Model** for a Line × Tester design under the original RCBD teaching assumption. It reads the simulated CSV from `data/` and walks through:
    
    *   *Exploratory Data Analysis:* Visualizing yields across treatments.

    *   *Fixed-Effects Arithmetic (Griffing/Line-Tester):* Step-by-step arithmetic to estimate combining abilities from cell averages.
    
    *   *Mixed-Effects Modeling:* Fitting the mixed model using `sommer` in R, extracting random GCA/SCA BLUPs and variance components under Low P stress (`env2`).
    
    *   *Breeding Values ($BV$) and Genotypic Values ($GV$):* Calculating parental additive contributions and hybrid commercial potentials.
    
    *   *Heritability Estimation:* Estimating narrow-sense, broad-sense, and family-mean heritabilities.

*   `multilocation_alpha_lattice.qmd` & `docs/multilocation_alpha_lattice.html`: New notebook for Stelio's real data. It explains why `multilocation.csv` fits an alpha-lattice layout, diagnoses the duplicated `GEN = 49` rows, reconstructs candidate incomplete-block labels from row order for the clean `NHACOONGO` location, and illustrates the alpha lattice using `aplot`.

***

## Required R Packages

To execute the `.qmd` notebooks, ensure the following R packages are installed:

```R
install.packages(c("tidyverse", "sommer", "knitr", "kableExtra", "aplot"))
```

*   `tidyverse` is utilized for data manipulation, reshaping, and exploratory plotting.

*   `sommer` is the standard quantitative genetics mixed-model package in R, running the mmer solver to partition variance components and calculate BLUPs.

*   `knitr` and `kableExtra` are used to render styled tables.

*   `aplot` is used to combine the alpha-lattice layout plot with a supporting replication-count panel.

***

## Compilation Instructions

This repository is configured as a Quarto project (see `_quarto.yml`), with `docs/` as the output directory and the notebooks registered for rendering. The `docs/` folder is also what GitHub Pages serves, so a fresh render publishes the site. To compile every notebook into `docs/` in one shot, run:

```bash
quarto render
```

To render a single notebook, pass its filename — the project-level `output-dir: docs` still applies:

```bash
quarto render gca_calculation.qmd
```

To render only the new alpha-lattice notebook:

```bash
quarto render multilocation_alpha_lattice.qmd
```

*(Alternatively, you can open the `.qmd` files in RStudio or Positron and run the code chunks interactively).*

***

## References

Bibliographic metadata for both works lives in [`references.bib`](references.bib) and is rendered automatically by Quarto in `gca_calculation.qmd` (see §1.3, Phenotype Model).

*   Comstock, R. E., & Robinson, H. F. (1948). The Components of Genetic Variance in Populations of Biparental Progenies and Their Use in Estimating the Average Degree of Dominance. *Biometrics*, 4(4), 254–266. [doi:10.2307/3001412](https://doi.org/10.2307/3001412)

*   Isik, F., Holland, J., & Maltecca, C. (2017). *Genetic Data Analysis for Plant and Animal Breeding*. Springer International Publishing. [doi:10.1007/978-3-319-55177-7](https://doi.org/10.1007/978-3-319-55177-7)
