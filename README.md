# Maize Combining Ability Analysis (GCA and SCA)

A highly structured R-based repository designed to instruct plant breeding students in calculating General Combining Ability (GCA) and Specific Combining Ability (SCA) from a maize dataset. 

This project is tailored specifically to Stélio Boaventura Nuvunga's PhD research on **"Genetic Analysis of Tolerance to Low Phosphorus in Maize"**.

***

## Breeding and Experimental Design

The primary genetic analysis focuses on Stelio's actual experimental mating and field layout:

1.  **Breeding Mating Scheme:** **7 × 7 Factorial Cross (Line × Tester / North Carolina Design II)**.

    *   **Females (Lines):** 7 Low-Phosphorus tolerant inbred donor lines sourced from CIMMYT (`CML364`, `CML366`, `CML434`, `CML435`, `CML439`, `CML440`, `CML532`).

    *   **Males (Testers):** 7 elite local Mozambican inbred lines (`MOZ1` to `MOZ7`) adapted to local agroecological zones.

    *   **Evaluated Progeny:** $F_1$ diallelic factorial hybrids ($N = 49$ unique families). (No selfs or reciprocals are represented in the primary analysis).

2.  **Field Trial Layout:** **Randomized Complete Block Design (RCBD)** with **3 replications (blocks)** nested within each environment.

3.  **Environments:** Combined evaluation across **2 environments** representing phosphorus treatments:

    *   `env1`: Normal/High Phosphorus (control).

    *   `env2`: Low Phosphorus Stress (treatment).

4.  **Evaluated Trait:** Grain Yield (tons/ha).

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
├── data/
│   └── maize_factorial_yield.csv      # Exported simulated dataset (294 plot observations)
└── docs/
    ├── iris_anova_cld.html            # Compiled HTML report for Iris ANOVA reference
    ├── diallel_gca_sca_analysis.html  # Compiled HTML report for Diallel GCA/SCA reference
    ├── gca_simulation.html            # Compiled HTML report for Maize Factorial simulation
    └── gca_calculation.html           # Compiled HTML report for Maize Factorial GCA/SCA calculations
```

### File Details

*   `data/maize_factorial_yield.csv`: Simulated grain yields containing 294 rows with columns: `env`, `block`, `female`, `male`, `cross`, and `yield`. The dataset incorporates inbreeding depression on the parental lines, deterministic GCA effects, specific heterosis (SCA), and GCA × Environment interactions (amplified CIMMYT tolerance under low phosphorus stress).

*   `gca_simulation.qmd` & `docs/gca_simulation.html`: The background script that sets up the parameters, builds the linear simulation model, and exports the data. Kept separate to keep Stelio's focus on calculations.

*   `gca_calculation.qmd` & `docs/gca_calculation.html`: Stelio's main learning notebook. It outlines the **Breeding Scheme**, **Experimental Design**, and **Phenotype Model** of his exact Line × Tester design (no complex diallel overlays needed). It reads the CSV from `data/` and walks through:
    
    *   *Exploratory Data Analysis:* Visualizing yields across treatments.

    *   *Fixed-Effects Arithmetic (Griffing/Line-Tester):* Step-by-step arithmetic to estimate combining abilities from cell averages.
    
    *   *Mixed-Effects Modeling:* Fitting the mixed model using `sommer` in R, extracting random GCA/SCA BLUPs and variance components under Low P stress (`env2`).
    
    *   *Breeding Values ($BV$) and Genotypic Values ($GV$):* Calculating parental additive contributions and hybrid commercial potentials.
    
    *   *Heritability Estimation:* Estimating narrow-sense, broad-sense, and family-mean heritabilities.

***

## Required R Packages

To execute the `.qmd` notebooks, ensure the following R packages are installed:

```R
install.packages(c("tidyverse", "sommer", "knitr", "kableExtra"))
```

*   `tidyverse` is utilized for data manipulation, reshaping, and exploratory plotting.

*   `sommer` is the standard quantitative genetics mixed-model package in R, running the mmer solver to partition variance components and calculate BLUPs.

*   `knitr` and `kableExtra` are used to render styled tables.

***

## Compilation Instructions

This repository is configured as a Quarto project (see `_quarto.yml`), with `docs/` as the output directory and all four notebooks registered for rendering. The `docs/` folder is also what GitHub Pages serves, so a fresh render publishes the site. To compile every notebook into `docs/` in one shot, run:

```bash
quarto render
```

To render a single notebook, pass its filename — the project-level `output-dir: docs` still applies:

```bash
quarto render gca_calculation.qmd
```

*(Alternatively, you can open the `.qmd` files in RStudio or Positron and run the code chunks interactively).*

***

## References

Bibliographic metadata for both works lives in [`references.bib`](references.bib) and is rendered automatically by Quarto in `gca_calculation.qmd` (see §1.3, Phenotype Model).

*   Comstock, R. E., & Robinson, H. F. (1948). The Components of Genetic Variance in Populations of Biparental Progenies and Their Use in Estimating the Average Degree of Dominance. *Biometrics*, 4(4), 254–266. [doi:10.2307/3001412](https://doi.org/10.2307/3001412)

*   Isik, F., Holland, J., & Maltecca, C. (2017). *Genetic Data Analysis for Plant and Animal Breeding*. Springer International Publishing. [doi:10.1007/978-3-319-55177-7](https://doi.org/10.1007/978-3-319-55177-7)
