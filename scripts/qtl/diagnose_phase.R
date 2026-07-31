#!/usr/bin/env Rscript
# qtl_diagnose_phase.R -- why did est.map return 161,170 cM instead of ~1,400?
#
# MY ERROR: I claimed recombination fractions are "orientation-invariant", so
# markers without a CML530 allele call could carry an arbitrary A/B label. That
# is false. Flipping BOTH markers of a pair leaves rf unchanged; flipping ONE
# turns rf into 1-rf. Per-marker polarity must be consistent IN LINKAGE PHASE
# along the chromosome or every adjacent pair looks unlinked -- which is exactly
# a ~100x inflated map.
#
# This script does NOT fix anything. It distinguishes three hypotheses before any
# rebuild, because the remedy differs completely:
#
#   H1 PHASE ONLY. Linkage is present but polarity is inconsistent.
#      Signature: min(rf, 1-rf) SMALL for physically adjacent markers, while rf
#      itself is scattered around 0.5. Remedy: infer phase from linkage.
#
#   H2 CML530 IS NOT THIS F2's PARENT (or the calls are uninformative).
#      Signature: markers polarized FROM CML530 are no better in phase than
#      markers polarized arbitrarily. If CML530 were a true parent, its calls
#      would put ~all of those markers in one consistent phase.
#      Remedy: stop using CML530 for polarity; ask Stelio which cross this is.
#
#   H3 SAMPLE IDENTITY IS BROKEN (the plant N <-> field-book risk).
#      Signature: min(rf, 1-rf) ALSO near 0.5 -- no linkage at any polarity.
#      Remedy: stop; get the plate submission sheet.
#
# Reads   data/qtl/derived/rqtl_gen_abh.csv, rqtl_phe.csv, abh_marker_info.tsv
# Writes  nothing.
# Usage: Rscript scripts/qtl/qtl_diagnose_phase.R > agent/qtl_diagnose_phase.log 2>&1

suppressMessages({library(qtl); library(data.table)})
rule <- function(x) cat("\n", strrep("=", 70), "\n", x, "\n", strrep("=", 70), "\n", sep = "")

D <- "data/qtl/derived"
cr <- read.cross(format = "csvsr", dir = "",
                 genfile = file.path(D, "rqtl_gen_abh.csv"),
                 phefile = file.path(D, "rqtl_phe.csv"),
                 genotypes = c("A", "H", "B"), na.strings = c("-", "NA"),
                 crosstype = "f2", estimate.map = FALSE)
mkinfo <- fread(file.path(D, "abh_marker_info.tsv"),
                colClasses = list(character = "chr_v5"))

rule("1. RAW TWO-LOCUS COUNTS FOR PHYSICALLY ADJACENT MARKERS")
# Estimate rf directly from genotype pairs, WITHOUT R/qtl's [0,0.5] constraint,
# so a flipped marker shows up as rf near 1 rather than being folded to 0.5.
G <- pull.geno(cr)                          # 1=AA 2=AB 3=BB
mk_chr <- rep(names(cr$geno), nmar(cr))
# expected recombinant fraction from an F2 two-locus genotype pair, additively:
# treat AA=0, AB=1, BB=2 dosage and use 1 - |corr| logic via a simple estimator
dose <- G - 1                               # 0,1,2
res <- rbindlist(lapply(names(cr$geno), function(ch) {
  idx <- which(mk_chr == ch); if (length(idx) < 2) return(NULL)
  d <- dose[, idx, drop = FALSE]
  n <- length(idx)
  r <- sapply(seq_len(n - 1), function(i) {
    ok <- is.finite(d[, i]) & is.finite(d[, i + 1])
    if (sum(ok) < 30) return(NA_real_)
    cor(d[ok, i], d[ok, i + 1])
  })
  data.table(chr = ch, i = seq_len(n - 1), marker = colnames(d)[-n], cor_next = r)
}))
# For tightly linked markers in correct phase, cor -> +1; flipped phase -> -1;
# unlinked -> 0. This separates the three hypotheses cleanly.
cat("correlation between physically adjacent markers (additive dosage):\n")
print(round(summary(res$cor_next), 3))
cat("\ndistribution:\n")
print(table(cut(res$cor_next, breaks = c(-1, -0.7, -0.3, 0.3, 0.7, 1),
                labels = c("<-0.7 (flipped, tight)", "-0.7..-0.3 (flipped)",
                           "-0.3..0.3 (UNLINKED)", "0.3..0.7 (phase ok)",
                           ">0.7 (phase ok, tight)"))))
cat("\n|cor| -- linkage strength regardless of phase:\n")
print(round(summary(abs(res$cor_next)), 3))

rule("2. VERDICT ON H3 (sample identity)")
frac_unlinked <- mean(abs(res$cor_next) < 0.3, na.rm = TRUE)
cat(sprintf("adjacent pairs with |cor| < 0.3 (no linkage at ANY phase): %.1f%%\n",
            100 * frac_unlinked))
if (frac_unlinked > 0.5) {
  cat("=> H3 LIKELY: linkage is largely absent. Phase cannot explain this.\n")
  cat("   STOP and verify the plant N <-> field-book correspondence.\n")
} else {
  cat("=> H3 REJECTED: strong linkage is present between physical neighbours.\n")
  cat("   The genotypes and the physical order are sound; only PHASE is wrong.\n")
}

rule("3. VERDICT ON H1 vs H2 (is CML530 polarity informative?)")
res <- merge(res, mkinfo[, .(marker, polarity_source, cml_class)],
             by = "marker", all.x = TRUE)
# A pair is "in phase" if cor > 0. Ask whether pairs whose FIRST marker was
# polarized from CML530 are more often in phase than arbitrary ones.
res[, first_from_cml := polarity_source == "CML530"]
tb <- res[!is.na(cor_next), .(n = .N,
                              pct_in_phase = round(100 * mean(cor_next > 0), 1),
                              median_cor = round(median(cor_next), 3)),
          by = first_from_cml]
print(tb)
cat("\nIf CML530 were a true parent of THIS F2, markers polarized from it would\n")
cat("be ~100%% in phase with each other. Arbitrary ones would be ~50%%.\n")

# sharper test: consecutive pairs where BOTH markers came from CML530
res[, both_cml := first_from_cml & shift(first_from_cml, -1, type = "lag") %in% TRUE]
bp <- res[!is.na(cor_next)]
bp[, prev_src := shift(polarity_source, 1), by = chr]
both <- bp[polarity_source == "CML530" & prev_src == "CML530"]
if (nrow(both) > 20) {
  cat(sprintf("\npairs with BOTH markers polarized from CML530: n=%d, %.1f%% in phase, median cor %.3f\n",
              nrow(both), 100 * mean(both$cor_next > 0), median(both$cor_next)))
}
arb <- bp[polarity_source == "arbitrary" | prev_src == "arbitrary"]
if (nrow(arb) > 20) {
  cat(sprintf("pairs involving an ARBITRARY marker      : n=%d, %.1f%% in phase, median cor %.3f\n",
              nrow(arb), 100 * mean(arb$cor_next > 0), median(arb$cor_next)))
}
cat("\nby CML530 class of the first marker:\n")
print(res[!is.na(cor_next), .(n = .N, pct_in_phase = round(100 * mean(cor_next > 0), 1)),
          by = cml_class][order(-n)])

rule("4. WHAT A CORRECT PHASE PASS WOULD RECOVER")
# Greedy phase fix along the physical order: flip marker i+1 whenever its
# correlation with the running reference is negative. This is the standard
# linkage-based phasing and needs NO parental genotype.
flips <- res[!is.na(cor_next), .(n_neg = sum(cor_next < 0), n = .N), by = chr]
flips[, pct_neg := round(100 * n_neg / n, 1)]
print(flips[order(as.integer(chr))])
cat("\nEach negative-correlation junction is one polarity flip to apply.\n")
cat("After phasing, |cor| is unchanged but all signs become positive, so rf\n")
cat("collapses to its true value and map length should fall to the maize range.\n")
cat("\nNOTE: phase can be inferred from LINKAGE ALONE. A parental genotype is\n")
cat("needed only to LABEL which resulting haplotype belongs to which parent --\n")
cat("i.e. for the SIGN of allele effects, not for map construction. I had this\n")
cat("backwards: I used CML530 to SET the coding, when it should only NAME it.\n")
