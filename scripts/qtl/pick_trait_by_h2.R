#!/usr/bin/env Rscript
# qtl_pick_trait_by_h2.R -- choose which F2 trait to map first, on evidence.
#
# THE PREMISE PROBLEM: heritability is NOT ESTIMABLE in this F2. One plant per
# genotype, zero within-genotype replication => genetic and micro-environmental
# variance cannot be separated. Any h2 computed from mozpue_phenotype.xlsx would
# be fabricated. This script therefore does NOT compute F2 heritability.
#
# WHAT IT DOES INSTEAD: takes h2 from the DIALLEL (data/multilocation.csv, 3 reps
# x 5 environments), which is replicated and where h2 IS estimable, and uses it
# as a PROXY ranking for which F2 traits are most likely to yield detectable QTL.
# Uses the identical sommer::mmer specification as variance_components.qmd so the
# numbers reconcile with the published table rather than being a second opinion.
#
# WHY A PROXY IS ONLY A PROXY -- the two experiments differ in:
#   * material  : F1 hybrids (diallel) vs F2 segregants
#   * scale     : per-plot over ~25 plants vs per-plant
#   * P regime  : 10 / 50 kg P/ha vs 20 kg P/ha
#   * genetics  : among-cross variance in a factorial vs within-cross segregation
# A trait can be heritable among hybrids and still segregate poorly in one F2,
# and vice versa. This ranking narrows the choice; it does not determine it.
#
# Reads   data/multilocation.csv via R/prep_traits.R, data/phenotype_dictionary.csv
# Writes  nothing (report only)
# Usage: Rscript scripts/qtl/qtl_pick_trait_by_h2.R > agent/qtl_pick_trait_by_h2.log 2>&1

suppressMessages({library(data.table)})
rule <- function(x) cat("\n", strrep("=", 70), "\n", x, "\n", strrep("=", 70), "\n", sep = "")

rule("0. PREFLIGHT")
need <- c("sommer", "tidyverse", "lme4")
inst <- rownames(installed.packages())
for (p in need) cat(sprintf("  %-12s %s\n", p, ifelse(p %in% inst, "OK", "MISSING")))
if (!"sommer" %in% inst)
  stop("sommer is required (variance_components.qmd uses sommer::mmer). ",
       "Install it, or this script cannot reproduce the published h2 values.")
suppressMessages({library(sommer); library(tidyverse)})

rule("1. WHICH TRAITS EXIST IN BOTH EXPERIMENTS")
dict <- fread("data/phenotype_dictionary.csv")
f2_traits    <- dict[project == "qtl_f2"        & role == "trait", variable]
multi_traits <- dict[project == "multilocation" & role == "trait", variable]
cat("F2 traits           :", paste(sort(f2_traits), collapse = " "), "\n")
cat("diallel traits      :", paste(sort(multi_traits), collapse = " "), "\n")
both <- sort(intersect(f2_traits, multi_traits))
cat("\nIN BOTH (mappable + h2 available):", paste(both, collapse = " "), "\n")
f2_only <- sort(setdiff(f2_traits, multi_traits))
cat("F2 ONLY (no h2 proxy possible)   :", paste(f2_only, collapse = " "), "\n")
m_only  <- sort(setdiff(multi_traits, f2_traits))
cat("DIALLEL ONLY (NOT measured in F2):", paste(m_only, collapse = " "), "\n")
cat("\n*** The diallel's most heritable traits are the FLOWERING traits\n")
cat("    (manuscript: DTA h2=0.60, DTS h2=0.61, Baker's ratio ~0.98).\n")
cat("    They are NOT in the F2 phenotype file, so they cannot be mapped here. ***\n")

rule("2. FIT THE DIALLEL MODEL (same spec as variance_components.qmd)")
source("R/prep_traits.R")
e_env <- 5; r_rep <- 3
cat("hybrids frame:", nrow(hybrids), "rows\n")
cat("analysis_traits:", paste(analysis_traits, collapse = " "), "\n\n")

fit_one <- function(tn, d) {
  d <- d %>% mutate(y = .data[[tn]]) %>% filter(!is.na(y))
  fit <- try(suppressMessages(suppressWarnings(mmer(
    y ~ env,
    random = ~ female + male + cross + female:env + male:env + cross:env +
               env:rep + env:rep:row + env:rep:col,
    data = d, verbose = FALSE))), silent = TRUE)
  if (inherits(fit, "try-error")) {
    cat(sprintf("  %-6s FIT FAILED\n", tn)); return(NULL)
  }
  vc <- summary(fit)$varcomp
  g <- function(k) { v <- vc$VarComp[sub("\\.y-y$", "", rownames(vc)) == k]
                     if (length(v)) max(0, v) else 0 }
  sL <- g("female"); sT <- g("male"); sS <- g("cross")
  sGE <- g("female:env") + g("male:env") + g("cross:env"); sE <- g("units")
  sG <- sL + sT + sS
  Vp <- sG + sGE / e_env + sE / (e_env * r_rep)
  out <- data.table(trait = tn, H2 = sG / Vp, h2 = (sL + sT) / Vp,
                    BR = if (2*(sL+sT) + sS > 0) 2*(sL+sT)/(2*(sL+sT)+sS) else NA_real_,
                    n = nrow(d))
  cat(sprintf("  %-6s H2=%.3f  h2=%.3f  BR=%s  n=%d\n", tn, out$H2, out$h2,
              ifelse(is.na(out$BR), "NA", sprintf("%.3f", out$BR)), out$n))
  out
}
res <- rbindlist(lapply(analysis_traits, fit_one, d = hybrids))

rule("3. RANKING, RESTRICTED TO F2-MAPPABLE TRAITS")
res[, mappable_in_f2 := trait %in% f2_traits]
setorder(res, -h2)
cat("ALL diallel traits by narrow-sense h2:\n")
print(res[, .(trait, H2 = round(H2, 3), h2 = round(h2, 3),
              BR = round(BR, 3), mappable_in_f2)])
cat("\nF2-MAPPABLE ONLY, by h2:\n")
print(res[mappable_in_f2 == TRUE, .(trait, H2 = round(H2, 3), h2 = round(h2, 3),
                                     BR = round(BR, 3))])

rule("4. RECOMMENDATION")
cand <- res[mappable_in_f2 == TRUE][order(-h2)]
top <- cand$trait[1]
cat(sprintf("Highest-h2 F2-mappable trait (diallel proxy): %s (h2=%.3f, H2=%.3f, BR=%.3f)\n",
            top, cand$h2[1], cand$H2[1], cand$BR[1]))
cat("\nCAVEATS THAT MUST TRAVEL WITH THIS CHOICE:\n")
cat("1. This is a PROXY from a different population, scale and P rate. It ranks\n")
cat("   plausibility, it does not predict F2 QTL detectability.\n")
cat("2. n=166 gives ~80% power only for QTL explaining >=8-10% of variance.\n")
cat("   Effects will be Beavis-inflated.\n")
derived <- dict[project == "qtl_f2" & role == "trait" & formula != "", variable]
cat("3. Derived/ratio traits in the F2 (", paste(derived, collapse = " "), ")\n", sep = "")
cat("   are FUNCTIONS of the measured traits and co-localize by arithmetic.\n")
if (top %in% derived)
  cat("   *** ", top, " IS a derived trait -- prefer mapping its components too. ***\n", sep = "")
cat("4. HI in the F2 is GYPP/SW (grain-per-shoot), NOT conventional harvest index,\n")
cat("   and the diallel's HI could not be reproduced from its own columns. An HI\n")
cat("   comparison across the two experiments is NOT like-for-like.\n")
cat("\nNOTE: no genetic map exists yet. est.map on the fixed B73 physical order,\n")
cat("with the distortion + find_quirky two-round pass, must run before any scan.\n")
