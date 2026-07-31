#!/usr/bin/env Rscript
# qtl_audits.R -- the ad-hoc audits I ran as inline `Rscript -e` blocks and then
# quoted numbers from. Written to a file so those numbers are reproducible.
# READ-ONLY.
# Usage: Rscript scripts/qtl/qtl_audits.R > agent/qtl_audits.log 2>&1

suppressMessages({library(qtl); library(data.table)})
rule <- function(x) cat("\n", strrep("=", 70), "\n", x, "\n", strrep("=", 70), "\n", sep = "")
D <- "data/qtl/derived"

rule("A. MAP POSITIONS -- the dummy 5 cM offsets that broke the first plot")
# est.map PRESERVES each chromosome's starting offset, and the input was
# read.cross's dummy 5 cM spacing. Minima came out as exact multiples of 5.
f <- file.path(D, "rqtl_cross_map.rds")
if (file.exists(f)) {
  m <- pull.map(readRDS(f))
  for (ch in names(m))
    cat(sprintf("chr %-3s n=%3d  min=%8.2f  max=%8.2f  span=%7.2f\n",
                ch, length(m[[ch]]), min(m[[ch]]), max(m[[ch]]),
                max(m[[ch]]) - min(m[[ch]])))
} else cat("absent:", f, "\n")

rule("B. WHAT 'A' MEANS IN THE PHASED FILE (the 1294 / 16 / 16 split)")
f <- file.path(D, "phase_flips.tsv")
if (file.exists(f)) {
  fl <- fread(f)
  cat("markers in phased file:", nrow(fl), "\n\n")
  cat("flipped x polarity_source:\n")
  print(dcast(fl[, .N, by = .(polarity_source, flipped)],
              polarity_source ~ flipped, value.var = "N", fill = 0))
  cat("\nflipped x cml_class:\n")
  print(dcast(fl[, .N, by = .(cml_class, flipped)],
              cml_class ~ flipped, value.var = "N", fill = 0))
  cat("\nA = CML530 allele                        :",
      fl[polarity_source == "CML530" & flipped == FALSE, .N], "\n")
  cat("A = TESTER allele (CML530 call overridden):",
      fl[polarity_source == "CML530" & flipped == TRUE, .N], "\n")
  cat("A = arbitrary, no CML530 call            :",
      fl[polarity_source == "arbitrary", .N], "\n")
} else cat("absent:", f, "\n")

rule("C. MARKERS WITH A = CML530 DEFINED (the 3,728)")
f <- file.path(D, "CML530_marker_alleles_direct.tsv")
if (file.exists(f)) {
  c530 <- fread(f)
  print(c530[, .N, by = cml_class][order(-N)])
  cat("\nA = CML530 DEFINED for:", c530[cml_class %in% c("ref","alt"), .N], "markers\n")
  cat("undefined (no-hit / third-allele):",
      c530[!cml_class %in% c("ref","alt"), .N], "\n")
} else cat("absent:", f, "\n")

rule("D. qc_pass BROKEN DOWN BY CRITERION -- how my 2,432-marker cut was composed")
d <- fread("data/qtl/Report_DMz26-3123_SNP_mapping_2.csv", skip = 6,
           header = TRUE, na.strings = c("", "NA"))
mc <- 1:which(names(d) == "RepAvg"); G <- as.matrix(d[, -mc, with = FALSE])
chr <- sub(" .*", "", d$`Chrom_Maize_B73_V4.0.assembly`)
n0 <- rowSums(G == "0"); n1 <- rowSums(G == "1"); n2 <- rowSums(G == "2")
tot <- n0 + n1 + n2
pref <- (2*n0 + n2) / (2*tot); maf <- pmin(pref, 1 - pref)
cs <- (n0 - tot/4)^2/(tot/4) + (n1 - tot/4)^2/(tot/4) + (n2 - tot/2)^2/(tot/2)
p121 <- pchisq(cs, 2, lower.tail = FALSE)
anch <- chr %in% as.character(1:10)
cat("anchored to chr1-10:", sum(anch), "of", nrow(d), "\n\n")
a <- data.table(maf = maf[anch], cr = d$CallRate[anch],
                rep = d$RepAvg[anch], p = p121[anch])
cat("--- each criterion applied ALONE to the anchored set ---\n")
cat(sprintf("MAF < 0.05  (monomorphic, unmappable) : %5d\n", a[maf < 0.05, .N]))
cat(sprintf("MAF < 0.15                            : %5d\n", a[maf < 0.15, .N]))
cat(sprintf("CallRate < 0.90                       : %5d\n", a[cr < 0.90, .N]))
cat(sprintf("RepAvg   < 0.95                       : %5d\n", a[rep < 0.95, .N]))
cat(sprintf("1:2:1 chi-sq p <= 0.01                : %5d\n", a[p <= 0.01, .N]))
b <- a[maf >= 0.15]
cat("\n--- among those with MAF >= 0.15 ---\n")
cat(sprintf("n                               : %5d\n", nrow(b)))
cat(sprintf("  also CallRate < 0.90          : %5d\n", b[cr < 0.90, .N]))
cat(sprintf("  also RepAvg < 0.95            : %5d\n", b[rep < 0.95, .N]))
cat(sprintf("  also 1:2:1 p <= 0.01          : %5d\n", b[p <= 0.01, .N]))
cat(sprintf("  passing all                   : %5d\n",
            b[cr >= 0.90 & rep >= 0.95 & p > 0.01, .N]))
cat(sprintf("  passing all EXCEPT 1:2:1      : %5d\n", b[cr >= 0.90 & rep >= 0.95, .N]))
