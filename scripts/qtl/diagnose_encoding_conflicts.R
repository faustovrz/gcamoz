#!/usr/bin/env Rscript
# diagnose_encoding_conflicts.R -- WHY do markers disagree with the A = CML530 coding?
#
# With A = CML530 at every marker, phase is fixed by construction: CML530's haplotype
# is A genome-wide. Adjacent markers must then correlate POSITIVELY. ~476 junctions in
# the current set saturate at rf = 0.5 instead. Three possible causes, and they have
# DIFFERENT signatures:
#
#   (1) WRONG PARENTAL CALL -- paralogous tag alignment, or the published assembly
#       differing from the actual CML530 plant. Signature: STRONG NEGATIVE
#       correlation (|cor| high, sign wrong). The marker is informative but mislabelled.
#   (2) NOISY F2 GENOTYPES -- heavy missing data or poor reproducibility. Signature:
#       correlation NEAR ZERO. No linkage at any polarity.
#   (3) genuine assembly/order error -- would show as a coherent block, not singletons.
#
# I asserted (1) without measuring. This separates them, and cross-tabs against the
# CallRate / RepAvg filters I dropped when going from 1,384 to 2,501 markers -- if the
# offenders are concentrated in low-CallRate markers, the cause is (2) and those
# filters were doing real work.
#
# READ-ONLY.
# Usage: Rscript scripts/qtl/diagnose_encoding_conflicts.R > agent/diagnose_encoding_conflicts.log 2>&1

suppressMessages({library(qtl); library(data.table)})
rule <- function(x) cat("\n", strrep("=", 70), "\n", x, "\n", strrep("=", 70), "\n", sep = "")
D <- "data/qtl/derived"

cr <- read.cross(format = "csvsr", dir = "",
                 genfile = file.path(D, "rqtl_gen_abh_all.csv"),
                 phefile = file.path(D, "rqtl_phe.csv"),
                 genotypes = c("A","H","B"), na.strings = c("-","NA"),
                 crosstype = "f2", estimate.map = FALSE)
info <- fread(file.path(D, "abh_all_marker_info.tsv"),
              colClasses = list(character = "chr_v5"))
d <- fread("data/qtl/Report_DMz26-3123_SNP_mapping_2.csv", skip = 6, header = TRUE)
qual <- d[, .(marker = AlleleID, CallRate, RepAvg)]

G <- pull.geno(cr); chrv <- rep(names(cr$geno), nmar(cr)); dose <- G - 1L
cat("markers:", ncol(dose), " individuals:", nrow(dose), "\n")

rule("1. PER-MARKER: BEST AND WORST CORRELATION WITH NEIGHBOURS")
K <- 5L
best <- worst <- rep(NA_real_, ncol(dose))
for (ch in unique(chrv)) {
  idx <- which(chrv == ch); n <- length(idx); dd <- dose[, idx, drop = FALSE]
  for (k in seq_len(n)) {
    nb <- setdiff(max(1,k-K):min(n,k+K), k)
    cs <- sapply(nb, function(j) {
      a <- dd[,k]; b <- dd[,j]; ok <- is.finite(a) & is.finite(b)
      if (sum(ok) < 30) return(NA_real_); cor(a[ok], b[ok])
    })
    if (all(is.na(cs))) next
    best[idx[k]]  <- max(cs, na.rm = TRUE)      # signed
    worst[idx[k]] <- min(cs, na.rm = TRUE)
  }
}
bt <- data.table(marker = colnames(dose), chr = chrv, best = best, worst = worst)
bt[, absmax := pmax(abs(best), abs(worst))]
# CLASSIFY
bt[, cls := fifelse(is.na(absmax), "uncomputable",
             fifelse(absmax < 0.30, "NOISY (no linkage any polarity)",
              fifelse(best > 0.30 & abs(worst) <= best, "OK (positive, in phase)",
               fifelse(abs(worst) > 0.30 & abs(worst) > best,
                       "MISLABELLED (strong negative)", "ambiguous"))))]
print(bt[, .N, by = cls][order(-N)])
cat("\nsigned best-correlation distribution:\n")
print(round(quantile(bt$best, c(0,.05,.10,.25,.50,.75,1), na.rm = TRUE), 3))

rule("2. CROSS-TAB AGAINST THE FILTERS I DROPPED")
bt <- merge(bt, qual, by = "marker", all.x = TRUE)
bt <- merge(bt, info[, .(marker, maf, p121)], by = "marker", all.x = TRUE)
print(bt[, .(n = .N,
             median_CallRate = round(median(CallRate, na.rm=TRUE), 3),
             median_RepAvg   = round(median(RepAvg, na.rm=TRUE), 3),
             median_maf      = round(median(maf, na.rm=TRUE), 3),
             pct_CallRate_lt_0.90 = round(100*mean(CallRate < 0.90, na.rm=TRUE), 1),
             pct_RepAvg_lt_0.95   = round(100*mean(RepAvg < 0.95, na.rm=TRUE), 1)),
         by = cls][order(-n)])

rule("3. WOULD CallRate / RepAvg HAVE REMOVED THE OFFENDERS?")
bad <- bt[cls %in% c("NOISY (no linkage any polarity)", "MISLABELLED (strong negative)")]
cat("offenders:", nrow(bad), sprintf("(%.1f%% of markers)\n", 100*nrow(bad)/nrow(bt)))
for (thr in c(0.80, 0.85, 0.90, 0.95)) {
  cat(sprintf("CallRate >= %.2f would remove %4d of %4d offenders and %4d of %4d sound markers\n",
              thr, bad[CallRate < thr, .N], nrow(bad),
              bt[!marker %in% bad$marker & CallRate < thr, .N],
              bt[!marker %in% bad$marker, .N]))
}
for (thr in c(0.90, 0.95, 0.98)) {
  cat(sprintf("RepAvg   >= %.2f would remove %4d of %4d offenders and %4d of %4d sound markers\n",
              thr, bad[RepAvg < thr, .N], nrow(bad),
              bt[!marker %in% bad$marker & RepAvg < thr, .N],
              bt[!marker %in% bad$marker, .N]))
}

rule("4. ARE THE MISLABELLED ONES CLUSTERED OR SCATTERED?")
mis <- bt[cls == "MISLABELLED (strong negative)"]
cat("mislabelled:", nrow(mis), "\n")
if (nrow(mis)) {
  print(mis[, .N, by = chr][order(as.integer(chr))])
  cat("\nthese are candidates for a WRONG CML530 CALL (paralogous tag, or the\n")
  cat("assembly differing from the parent plant). Scattered => sporadic; clustered\n")
  cat("=> a structural or assembly problem in that region.\n")
}

rule("VERDICT")
nz <- bt[cls == "NOISY (no linkage any polarity)", .N]
nm <- nrow(mis)
cat(sprintf("NOISY (bad F2 genotypes)        : %4d\n", nz))
cat(sprintf("MISLABELLED (bad CML530 call)   : %4d\n", nm))
cat("\nIf NOISY dominates, the cause is data quality and CallRate/RepAvg were doing\n")
cat("real work -- I was wrong to call them arbitrary. If MISLABELLED dominates, the\n")
cat("CML530 allele read is the problem at those markers and they should be dropped.\n")
