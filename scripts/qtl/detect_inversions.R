#!/usr/bin/env Rscript
# qtl_detect_inversions.R -- find segments where the B73 physical order is the WRONG
# order for this population, i.e. segregating inversions.
#
# WHY (FVRZ): chr4 and chr9 carry the well-known inv4m and inv9f inversion
# polymorphisms. If one parent carries the inverted haplotype and the other the
# standard arrangement, then INSIDE the inverted segment the B73 marker order is
# reversed relative to this population. Building the map on the fixed B73 physical
# order therefore forces a reversed segment into forward order, and physically
# "adjacent" markers look highly recombinant.
#
# That is a STRUCTURAL explanation for what no amount of marker filtering fixed:
# chr4 stayed at 1,182 cM through every filter, while chr8 -- with MORE weak
# junctions -- sat at 125 cM. I kept treating it as marker quality. It is not.
#
# SIGNATURE: within a correctly ordered segment, |cor| decays with distance, so the
# pairwise matrix is diagonal-dominant. Within a segment whose order is REVERSED
# relative to the truth, the correlation structure runs ANTI-DIAGONAL: the first
# marker of the block correlates best with the last, and so on.
#
# TEST (order-agnostic, no literature coordinates assumed): slide a window along
# each chromosome and ask whether REVERSING the marker order inside it increases
# the sum of adjacent-pair correlations. A real inversion gives a large, contiguous
# window where reversal helps; noise gives small scattered gains.
#
# NOTE ON BIOLOGY: in an F2 segregating for an inversion, recombination is
# suppressed in inversion heterozygotes, so the segment also behaves as a single
# linkage block. Detecting the boundaries is the first step; how to HANDLE it
# (reverse the order, or treat the block as one locus) is a separate decision.
#
# READ-ONLY.
# Usage: Rscript scripts/qtl/qtl_detect_inversions.R > agent/qtl_detect_inversions.log 2>&1

suppressMessages({library(qtl); library(data.table)})
rule <- function(x) cat("\n", strrep("=", 70), "\n", x, "\n", strrep("=", 70), "\n", sep = "")

D <- "data/qtl/derived"
cr <- read.cross(format = "csvsr", dir = "",
                 genfile = file.path(D, "rqtl_gen_abh_phased.csv"),
                 phefile = file.path(D, "rqtl_phe.csv"),
                 genotypes = c("A","H","B"), na.strings = c("-","NA"),
                 crosstype = "f2", estimate.map = FALSE)
mk <- fread(file.path(D, "abh_marker_info.tsv"), colClasses = list(character = "chr_v5"))
G <- pull.geno(cr); chrv <- rep(names(cr$geno), nmar(cr)); dose <- G - 1L
cat("individuals:", nrow(dose), " markers:", ncol(dose), "\n")

corm <- function(d) {
  n <- ncol(d); m <- matrix(NA_real_, n, n)
  for (i in seq_len(n)) for (j in i:n) {
    a <- d[, i]; b <- d[, j]; ok <- is.finite(a) & is.finite(b)
    v <- if (sum(ok) < 30) NA_real_ else cor(a[ok], b[ok])
    m[i, j] <- v; m[j, i] <- v
  }
  m
}
adj_sum <- function(m) { n <- nrow(m); sum(abs(m[cbind(1:(n-1), 2:n)]), na.rm = TRUE) }

rule("1. SLIDING-WINDOW REVERSAL TEST")
# For each window, compare the sum of adjacent |cor| in the current order against
# the same sum with the window's markers reversed. Positive gain => the current
# order is worse than reversed => candidate inversion.
WINS <- c(10, 20, 40, 60, 80, 120)
hits <- rbindlist(lapply(names(cr$geno), function(ch) {
  idx <- which(chrv == ch); n <- length(idx)
  d <- dose[, idx, drop = FALSE]
  pos <- mk[match(colnames(d), marker), pos_v5]
  rbindlist(lapply(WINS, function(w) {
    if (n < w + 4) return(NULL)
    rbindlist(lapply(seq_len(n - w + 1), function(s) {
      e <- s + w - 1
      seg <- s:e
      # local sub-order: current vs reversed, evaluated on adjacent pairs only
      cur <- d[, seg, drop = FALSE]
      rev_ <- d[, rev(seg), drop = FALSE]
      g <- adj_sum(corm(rev_)) - adj_sum(corm(cur))
      data.table(chr = ch, win = w, start = s, end = e,
                 start_Mb = pos[s]/1e6, end_Mb = pos[e]/1e6, gain = g)
    }))
  }))
}))
cat("windows evaluated:", nrow(hits), "\n")
cat("\ngain distribution (reversed minus current, adjacent |cor| sum):\n")
print(round(summary(hits$gain), 3))
cat("\ntop 20 windows by gain:\n")
print(head(hits[order(-gain)], 20))

rule("2. PER-CHROMOSOME BEST CANDIDATE")
best <- hits[, .SD[which.max(gain)], by = chr]
setorder(best, -gain)
print(best[, .(chr, win, start_Mb = round(start_Mb,1), end_Mb = round(end_Mb,1),
               gain = round(gain, 2))])
cat("\nA real inversion should show a LARGE positive gain on chr4 and/or chr9 and\n")
cat("near-zero gains on the well-behaved chromosomes.\n")

rule("3. ANTI-DIAGONAL STRUCTURE ON chr4 AND chr9")
# Direct look at the signature: for each marker, is its best correlate a NEARBY
# marker (normal) or a DISTANT one (order wrong)?
for (ch in c("4", "9", "8")) {          # chr8 as a well-behaved control
  idx <- which(chrv == ch); d <- dose[, idx, drop = FALSE]
  pos <- mk[match(colnames(d), marker), pos_v5]
  m <- corm(d); diag(m) <- NA
  bestj <- apply(abs(m), 1, function(x) if (all(is.na(x))) NA_integer_ else which.max(x))
  lag <- abs(bestj - seq_along(bestj))
  cat(sprintf("\nchr%-2s n=%3d  |index lag to best correlate|: median %.0f, q90 %.0f, max %.0f\n",
              ch, length(idx), median(lag, na.rm = TRUE),
              quantile(lag, .9, na.rm = TRUE), max(lag, na.rm = TRUE)))
  far <- which(lag > 20)
  if (length(far)) {
    cat("  markers whose best correlate is >20 positions away:", length(far), "\n")
    cat("  their physical positions (Mb):",
        paste(round(pos[far]/1e6, 1), collapse = " "), "\n")
    cat("  positions of those best correlates (Mb):",
        paste(round(pos[bestj[far]]/1e6, 1), collapse = " "), "\n")
  } else cat("  none -- best correlate is always local (order looks right)\n")
}

rule("VERDICT")
cat("If chr4/chr9 show a block of markers whose best correlate is far away in\n")
cat("INDEX but nearby in the REVERSED index, that is the inversion. The fix is to\n")
cat("reverse the marker order within the segment before est.map -- not to delete\n")
cat("markers. Deleting them, which is what I have been doing, throws away real\n")
cat("data and cannot straighten the order.\n")
cat("\nCAUTION: reversal gain can also arise from a genuine assembly error in B73,\n")
cat("or from a translocation. Reversal is the correct fix for all of those; the\n")
cat("biological label (inv4m / inv9f) needs the literature coordinates to confirm.\n")
