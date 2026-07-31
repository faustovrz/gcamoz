#!/usr/bin/env Rscript
# verify_dart_codes.R -- prove which DArT one-row code is the heterozygote, directly.
#
# The DArT report ships its OWN summary columns: FreqHomRef, FreqHomSnp, FreqHets.
# The per-marker frequency of each genotype code must match the corresponding
# column exactly. That is definitive, unlike inferring it from the F2 expectation.
#
# Under test:
#   code "0" == hom reference   -> should match FreqHomRef
#   code "1" == hom SNP/alt     -> should match FreqHomSnp
#   code "2" == HETEROZYGOTE    -> should match FreqHets
# The dosage misreading would instead make "1" the het and "2" the hom-alt.
# READ-ONLY.
# Usage: Rscript scripts/qtl/verify_dart_codes.R > agent/verify_dart_codes.log 2>&1

suppressMessages(library(data.table))
rule <- function(x) cat("\n", strrep("=", 70), "\n", x, "\n", strrep("=", 70), "\n", sep = "")

d <- fread("data/qtl/Report_DMz26-3123_SNP_mapping_2.csv", skip = 6,
           header = TRUE, na.strings = c("", "NA"))
mc <- 1:which(names(d) == "RepAvg")
G  <- as.matrix(d[, -mc, with = FALSE])
cat("markers:", nrow(G), " samples:", ncol(G), "\n")

n0 <- rowSums(G == "0"); n1 <- rowSums(G == "1")
n2 <- rowSums(G == "2"); nm <- rowSums(G == "-")
called <- n0 + n1 + n2
# DArT's Freq* columns are fractions of CALLED genotypes
f0 <- n0 / called; f1 <- n1 / called; f2 <- n2 / called

cmp <- function(label, obs, col) {
  ref <- d[[col]]
  ok <- is.finite(obs) & is.finite(ref)
  mx <- max(abs(obs[ok] - ref[ok]))
  cat(sprintf("  %-28s vs %-12s  n=%d  max|diff| = %-12.3g  %s\n",
              label, col, sum(ok), mx, ifelse(mx < 1e-6, "MATCH", "NO MATCH")))
  mx < 1e-6
}

rule("1. THE CORRECT HYPOTHESIS: 0=homRef, 1=homSnp, 2=HET")
a <- cmp("freq of code 0", f0, "FreqHomRef")
b <- cmp("freq of code 1", f1, "FreqHomSnp")
c <- cmp("freq of code 2", f2, "FreqHets")

rule("2. THE DOSAGE MISREADING: 1=het, 2=homSnp")
cmp("freq of code 1", f1, "FreqHets")
cmp("freq of code 2", f2, "FreqHomSnp")

rule("3. CALL RATE CROSS-CHECK")
cmp("called / total samples", called / ncol(G), "CallRate")

rule("4. OneRatio COLUMNS (independent allele-frequency check)")
# OneRatioRef = frequency of the reference allele call among samples scoring it
cat("  columns present:", paste(intersect(c("OneRatioRef","OneRatioSnp"), names(d)),
                                collapse = ", "), "\n")

rule("VERDICT")
if (a && b && c) {
  cat("CONFIRMED by DArT's own summary columns:\n")
  cat("  code 0 = homozygous REFERENCE\n")
  cat("  code 1 = homozygous SNP (alternate)\n")
  cat("  code 2 = HETEROZYGOTE\n")
  cat("\nNOT a minor-allele dosage. A dosage reading would put the het at 1 and\n")
  cat("hom-alt at 2, systematically swapping H and B on every marker.\n")
  cat("So `ABH[Gm == \"2\"] <- \"H\"` is correct.\n")
} else {
  cat("NOT CONFIRMED -- the encoding assumption is WRONG. See section 2.\n")
}
