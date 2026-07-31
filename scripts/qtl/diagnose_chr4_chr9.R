#!/usr/bin/env Rscript
# qtl_diagnose_chr4_chr9.R -- localise what inflates chr4 (1,274 cM) and chr9 (407 cM)
# when the other eight chromosomes sit at a plausible 118-233 cM.
#
# DO NOT TUNE THRESHOLDS TO HIT A TARGET LENGTH. The point of this script is to
# identify WHICH markers break these two chromosomes and WHY, so the filter is
# chosen on evidence rather than fitted to the maize expectation.
#
# Large est.map gaps correspond to low-|cor| junctions, so the offenders can be
# found straight from the dosage matrix without re-running est.map.
#
# QUESTIONS
#   Q1 Are the offenders SINGLETONS or a contiguous BLOCK? A block in the wrong
#      place is a liftover/misassembly artefact (the chr7 v2->v5 case in TeoNAM);
#      scattered singletons are paralogs or genotyping error.
#   Q2 Do they coincide with a physical feature -- a huge Mb gap, a chromosome end,
#      the pericentromere?
#   Q3 What does the local-linkage filter see in them? If they sit just above
#      LINK_MIN = 0.30, the threshold is simply too permissive for a map this
#      dense (median |cor| 0.970, Q1 0.951).
#   Q4 Are they enriched for arbitrary polarity, no-hit, or low CallRate?
#
# READ-ONLY.
# Usage: Rscript scripts/qtl/qtl_diagnose_chr4_chr9.R > agent/qtl_diagnose_chr4_chr9.log 2>&1

suppressMessages({library(qtl); library(data.table)})
rule <- function(x) cat("\n", strrep("=", 70), "\n", x, "\n", strrep("=", 70), "\n", sep = "")

D <- "data/qtl/derived"
cr <- read.cross(format = "csvsr", dir = "",
                 genfile = file.path(D, "rqtl_gen_abh_phased.csv"),
                 phefile = file.path(D, "rqtl_phe.csv"),
                 genotypes = c("A","H","B"), na.strings = c("-","NA"),
                 crosstype = "f2", estimate.map = FALSE)
mk <- fread(file.path(D, "abh_marker_info.tsv"), colClasses = list(character = "chr_v5"))
fl <- fread(file.path(D, "phase_flips.tsv"), colClasses = list(character = "chr"))
G <- pull.geno(cr); chrv <- rep(names(cr$geno), nmar(cr)); dose <- G - 1L

# junction table: |cor| between physically adjacent markers, per chromosome
jt <- rbindlist(lapply(names(cr$geno), function(ch) {
  i <- which(chrv == ch); if (length(i) < 2) return(NULL)
  d <- dose[, i, drop = FALSE]; n <- length(i)
  r <- sapply(seq_len(n-1), function(k) {
    a <- d[,k]; b <- d[,k+1]; ok <- is.finite(a) & is.finite(b)
    if (sum(ok) < 30) return(NA_real_); cor(a[ok], b[ok])
  })
  data.table(chr = ch, k = seq_len(n-1),
             m1 = colnames(d)[-n], m2 = colnames(d)[-1], cor = r)
}))
jt <- merge(jt, mk[, .(m1 = marker, pos1 = pos_v5)], by = "m1")
jt <- merge(jt, mk[, .(m2 = marker, pos2 = pos_v5)], by = "m2")
jt[, mb_gap := (pos2 - pos1) / 1e6]
setorder(jt, chr, k)

rule("1. JUNCTION QUALITY PER CHROMOSOME")
print(jt[, .(n = .N,
             median_abs_cor = round(median(abs(cor), na.rm = TRUE), 3),
             n_weak_lt0.7 = sum(abs(cor) < 0.7, na.rm = TRUE),
             n_weak_lt0.3 = sum(abs(cor) < 0.3, na.rm = TRUE),
             n_negative = sum(cor < 0, na.rm = TRUE)),
         by = chr][order(as.integer(chr))])
cat("\n=> compare chr4 and chr9 against the eight well-behaved chromosomes.\n")

rule("2. Q1/Q2 -- WEAK JUNCTIONS ON chr4 AND chr9: SINGLETONS OR BLOCKS?")
for (ch in c("4", "9")) {
  cat("\n---------------- chr", ch, " ----------------\n", sep = "")
  w <- jt[chr == ch & abs(cor) < 0.7][order(k)]
  cat("weak junctions (|cor| < 0.7):", nrow(w), "of", jt[chr == ch, .N], "\n")
  if (nrow(w)) {
    print(w[, .(k, cor = round(cor, 3), mb_gap = round(mb_gap, 2),
                pos1_Mb = round(pos1/1e6, 1))][1:min(25, .N)])
    d <- diff(w$k)
    cat("consecutive-index runs (1 = contiguous block):\n")
    print(table(d)[1:min(6, length(table(d)))])
    cat("longest contiguous run of weak junctions:",
        max(rle(c(TRUE, d == 1))$lengths), "\n")
  }
}

rule("3. Q2 -- IS IT A PHYSICAL GAP OR A CHROMOSOME END?")
print(jt[, .(max_mb_gap = round(max(mb_gap, na.rm = TRUE), 1),
             median_mb_gap = round(median(mb_gap, na.rm = TRUE), 2)),
         by = chr][order(as.integer(chr))])
cat("\nweak junctions vs Mb gap size (all chromosomes):\n")
jt[, gapbin := cut(mb_gap, c(-Inf, 0.5, 2, 10, 30, Inf))]
print(jt[!is.na(cor), .(n = .N, pct_weak = round(100*mean(abs(cor) < 0.7), 1)), by = gapbin][order(gapbin)])

rule("4. Q3 -- WHERE DO THE OFFENDERS SIT RELATIVE TO LINK_MIN = 0.30?")
# per-marker max |cor| with neighbours, same computation as the filter
K <- 5L
best <- rep(NA_real_, ncol(dose))
for (ch in unique(chrv)) {
  idx <- which(chrv == ch); n <- length(idx); dd <- dose[, idx, drop = FALSE]
  for (k in seq_len(n)) {
    nb <- setdiff(max(1,k-K):min(n,k+K), k)
    cs <- sapply(nb, function(j) {
      a <- dd[,k]; b <- dd[,j]; ok <- is.finite(a) & is.finite(b)
      if (sum(ok) < 30) return(NA_real_); abs(cor(a[ok], b[ok]))
    })
    best[idx[k]] <- suppressWarnings(max(cs, na.rm = TRUE))
  }
}
bt <- data.table(marker = colnames(dose), chr = chrv, best = best)
cat("distribution of per-marker max |cor|:\n")
print(round(quantile(bt$best, c(0,.01,.02,.05,.10,.25,.50,.75,1), na.rm = TRUE), 3))
cat("\ncount below each candidate threshold, ALL chromosomes:\n")
for (t in c(0.3, 0.5, 0.6, 0.7, 0.8, 0.9))
  cat(sprintf("  |cor| < %.1f : %4d markers (%.1f%%)\n", t,
              sum(bt$best < t, na.rm = TRUE), 100*mean(bt$best < t, na.rm = TRUE)))
cat("\nsame, chr4 and chr9 only:\n")
for (t in c(0.3, 0.5, 0.6, 0.7, 0.8, 0.9))
  cat(sprintf("  |cor| < %.1f : chr4 %3d, chr9 %3d\n", t,
              bt[chr == "4" & best < t, .N], bt[chr == "9" & best < t, .N]))

rule("5. Q4 -- ARE OFFENDERS ENRICHED FOR ANY MARKER PROPERTY?")
bt <- merge(bt, mk[, .(marker, CallRate, maf, p121, cml_class, polarity_source)],
            by = "marker", all.x = TRUE)
bt <- merge(bt, fl[, .(marker, flipped)], by = "marker", all.x = TRUE)
bt[, weak := best < 0.7]
cat("by cml_class:\n")
print(bt[, .(n = .N, pct_weak = round(100*mean(weak, na.rm = TRUE), 1)), by = cml_class][order(-n)])
cat("\nby whether phasing flipped the marker:\n")
print(bt[, .(n = .N, pct_weak = round(100*mean(weak, na.rm = TRUE), 1)), by = flipped])
cat("\nCallRate / MAF, weak vs sound:\n")
print(bt[!is.na(weak), .(n = .N,
                         median_CallRate = round(median(CallRate), 3),
                         median_maf = round(median(maf), 3),
                         median_p121 = signif(median(p121), 3)), by = weak])

rule("VERDICT")
cat("Read section 2 first: a long contiguous run of weak junctions means a\n")
cat("DISPLACED BLOCK (wrong physical position) and the fix is to drop or relocate\n")
cat("that block. Scattered singletons mean paralogous or erroneous markers and the\n")
cat("fix is a stricter local-linkage threshold -- justified by section 4 only if\n")
cat("the offenders sit in a clear tail, not if raising the threshold would start\n")
cat("removing markers from the well-behaved chromosomes too.\n")
