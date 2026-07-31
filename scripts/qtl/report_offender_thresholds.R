#!/usr/bin/env Rscript
# report_offender_thresholds.R -- offender count across a range of cutoffs, so the
# threshold is FVRZ's choice and not mine.
#
# DEFINITION. A = CML530 at every marker, so phase is fixed by construction and a
# marker must correlate POSITIVELY with its physical neighbours. A marker that
# correlates with NONE of them contradicts the encoding, which means the marker is
# wrong -- not mis-phased. That principle is fixed; only the cutoff is a choice.
#
# Two knobs, both reported rather than assumed:
#   thr -- the |cor| a marker must reach with at least ONE neighbour
#   K   -- how many neighbours each side to consider
#
# Cannot be reached by any DArT statistic -- measured, not assumed:
#   RepAvg  offender rate is FLAT across bands (4.1% at the bottom vs 3.6% at the top)
#   CallRate has a real gradient but 80 of 117 offenders sit above 0.95
#   MAF, co-tag concordance: also no
# (agent/report_callrate_repavg.log, agent/check_tag_redundancy.log)
#
# READ-ONLY. Writes nothing, changes nothing.
# Usage: Rscript scripts/qtl/report_offender_thresholds.R > agent/report_offender_thresholds.log 2>&1

suppressMessages({library(qtl); library(data.table)})
rule <- function(x) cat("\n", strrep("=", 70), "\n", x, "\n", strrep("=", 70), "\n", sep = "")
D <- "data/qtl/derived"

cr <- read.cross(format = "csvsr", dir = "",
                 genfile = file.path(D, "rqtl_gen_abh_all.csv"),
                 phefile = file.path(D, "rqtl_phe.csv"),
                 genotypes = c("A","H","B"), na.strings = c("-","NA"),
                 crosstype = "f2", estimate.map = FALSE)
G <- pull.geno(cr); chrv <- rep(names(cr$geno), nmar(cr)); dose <- G - 1L
NM <- colnames(dose)
cat("markers:", length(NM), " individuals:", nrow(dose), "\n")

# per-marker max |cor| with neighbours, for each K
maxcor <- function(K) {
  out <- rep(NA_real_, ncol(dose))
  for (ch in unique(chrv)) {
    idx <- which(chrv == ch); n <- length(idx); dd <- dose[, idx, drop = FALSE]
    for (k in seq_len(n)) {
      nb <- setdiff(max(1, k - K):min(n, k + K), k)
      cs <- sapply(nb, function(j) {
        a <- dd[, k]; b <- dd[, j]; ok <- is.finite(a) & is.finite(b)
        if (sum(ok) < 30) return(NA_real_)
        abs(cor(a[ok], b[ok]))
      })
      if (!all(is.na(cs))) out[idx[k]] <- max(cs, na.rm = TRUE)
    }
  }
  out
}
KS <- c(1L, 2L, 3L, 5L, 10L)
MC <- lapply(KS, maxcor); names(MC) <- paste0("K", KS)

rule("1. DISTRIBUTION OF max |cor| WITH ANY NEIGHBOUR")
for (i in seq_along(KS)) {
  cat(sprintf("\nK = %-2d  ", KS[i]))
  print(round(quantile(MC[[i]], c(0,.01,.02,.03,.05,.10,.25,.50,1), na.rm = TRUE), 3))
}
cat("\nThe cliff between the 2nd and 5th percentile is where the offenders separate\n")
cat("from the bulk. A cutoff inside that gap is insensitive; one above it starts\n")
cat("taking sound markers.\n")

rule("2. OFFENDER COUNT ACROSS THRESHOLDS  (rows = |cor| cutoff, cols = K)")
THR <- c(0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.50, 0.60, 0.70, 0.80)
tab <- data.table(threshold = THR)
for (i in seq_along(KS))
  tab[[paste0("K", KS[i])]] <- sapply(THR, function(t) sum(MC[[i]] < t, na.rm = TRUE))
print(tab)
cat("\nsame, as % of", length(NM), "markers:\n")
pct <- copy(tab)
for (nmc in setdiff(names(pct), "threshold"))
  pct[[nmc]] <- round(100 * pct[[nmc]] / length(NM), 2)
print(pct)

rule("3. WHERE IS THE COUNT INSENSITIVE TO THE CUTOFF?")
cat("marginal markers added per 0.05 step of threshold (K = 5):\n")
m5 <- MC[[which(KS == 5L)]]
st <- seq(0.10, 0.80, by = 0.05)
mg <- data.table(from = head(st, -1), to = tail(st, -1))
mg[, added := sapply(seq_len(.N), function(i)
  sum(m5 < to[i], na.rm = TRUE) - sum(m5 < from[i], na.rm = TRUE))]
print(mg)
cat("\nA plateau of near-zero `added` means the threshold is in the gap and the\n")
cat("choice does not matter. Where `added` starts climbing, the cutoff has reached\n")
cat("the bulk of sound markers.\n")

rule("4. PER-CHROMOSOME AT A FEW CUTOFFS (K = 5)")
for (t in c(0.20, 0.30, 0.40)) {
  cat(sprintf("\n|cor| < %.2f  -> %d markers\n", t, sum(m5 < t, na.rm = TRUE)))
  print(table(chrv[which(m5 < t)]))
}

rule("5. OVERLAP WITH THE ROUTE-DISAGREEMENT MARKERS")
f <- file.path(D, "CML530_marker_alleles_direct.tsv")
cat("The 36 markers where the two coordinate routes disagreed are flagged\n")
cat("independently (verify_route_convergence.R). Whether they overlap the offenders\n")
cat("is worth knowing -- if they do, one filter handles both.\n")
cat("(not computed here: the convergence output is not persisted to a file)\n")

rule("SUMMARY")
cat(sprintf("At K = 5: |cor| < 0.20 -> %d,  < 0.30 -> %d,  < 0.40 -> %d markers\n",
            sum(m5 < 0.20, na.rm = TRUE), sum(m5 < 0.30, na.rm = TRUE),
            sum(m5 < 0.40, na.rm = TRUE)))
cat("Your call on the cutoff and on K. Nothing has been removed.\n")
