#!/usr/bin/env Rscript
# qtl_qc_teonam.R -- STEP 2: FVRZ's QC, in his order, on the fully-encoded set.
#
# Input: data/qtl/derived/rqtl_gen_abh_all.csv -- 3,728 markers, 186 individuals,
# A = CML530 on every marker, NO filters applied (qtl_encode_abh_all.R).
#
# His pipeline, verified from zealhmm/scripts/teonam_qtl_permap.R:
#   1. segregation distortion: chi2 vs cross-type expected freqs, then
#      PER-CHROMOSOME renorm_z outliers, z > 1.96                       (lines 134-147)
#   1b. LD PRUNE at Pearson r = 0.95 (replaces findDupMarkers, which found
#      only 28 markers). See the block below for the derivation.
#   2. ROUND 1 est.map, error.prob = ERROR_PROB, map.function = "haldane" (lines 95-98, 150)
#      (the source pipeline used 0.001, its BC1S4 value; ERROR_PROB is set below)
#   3. find_quirky on the round-1 map: fine_thr = smallest outlier gap,
#      island_thr = 2 cM, island_max_n = 20                            (lines 158-170)
#   4. ROUND 2 est.map                                                 (line 177)
#   5. guardrails: drop rate > 10%; length 1348-1596 cM                (lines 188-201)
#
# Only adaptation: expected genotype frequencies BC1S4 (0.734, 0.031, 0.234) -> F2
# (0.25, 0.50, 0.25). Nothing else changed, nothing of mine added.
#
# Usage: Rscript scripts/qtl/qtl_qc_teonam.R > agent/qtl_qc_teonam.log 2>&1

suppressMessages({library(qtl); library(data.table); library(ggplot2)})
rule <- function(x) cat("\n", strrep("=", 70), "\n", x, "\n", strrep("=", 70), "\n", sep = "")

D <- "data/qtl/derived"
MAPTOOLS <- "/Users/fvrodriguez/repos/zealhmm/scripts/map_tools.R"
# 1% genotype error, set by FVRZ 2026-07-31 and used THROUGHOUT -- map construction here
# and calc.genoprob at scan time. Re-swept on the current marker set after the CML530
# reorder, the MAPQ/XA filter and the r = 0.95 prune (agent/sweep_error_prob.log):
#
#   error.prob   total cM   median gap   expected   ratio   max gap   cov
#     2.5e-3       3397.2     2.329        1.295     1.80     30.0    96.7%
#     5e-3         2445.8     1.195        1.292     0.93     29.7    96.7%
#     1e-2         1858.4     0.834        1.289     0.65     29.0    96.7%
#     2e-2         1552.0     0.658        1.286     0.51     27.8    97.9%
#
# The two criteria still disagree: median-gap calibration points at 5e-3 (ratio 0.93),
# the 1348-1596 target band at 2e-2. 1e-2 sits between them. The 1% figure is also what
# the map-inflation arithmetic implies independently (HANDOVER section 1.8).
#
# WHAT THE RE-SWEEP CHANGED versus the original one: marker count is now stable across
# the range (1,145-1,163) instead of moving with the parameter, coverage IMPROVES with
# error.prob (96.7 -> 98.0%) instead of collapsing, and max_gap now responds (30.3 ->
# 26.8) where it was previously immune. So the old warning that high values "crush the
# bulk by deleting markers" no longer applies -- only the median-gap ratio does.
#
# AT SCAN TIME it makes no measurable difference: PUE scanned at 0.01 vs 2.5e-3 on this
# map gave r = 0.9978 over 3,005 positions, max |dLOD| 0.30, mean +0.009
# (agent/scan_pue_errorprob.log). Use 0.01 in calc.genoprob for consistency.
ERROR_PROB <- 0.01; MAP_FUN <- "haldane"
# LD prune threshold, chosen by FVRZ 2026-07-31 in Pearson r -- the natural unit of the
# Teng & Xu estimator. Markers are kept pairwise >= -50*ln(r) cM apart on the
# correlation axis. Replaces findDupMarkers (see block 1b).
PRUNE_R <- 0.95; THIN_CM <- -50 * log(PRUNE_R)
# find_quirky parameters set by FVRZ for THIS map's scale. The TeoNAM defaults were
# invalid here, in opposite directions:
#   island_thr = 2 cM  -- BELOW our median gap of 2.308, so `diff(pos) > thr` was TRUE
#                         for over half of all junctions. That fragmented each
#                         chromosome into size-1/2 components, every one of which
#                         satisfied "isolated and small" -> 78.9% of markers flagged.
#                         TeoNAM's own comment puts 2 cM ABOVE its 99.99th-percentile
#                         gap (~1 cM); here it is below the median, off by ~50-100x.
#   fine_thr (data-driven, smallest outlier gap) -- came out 75.734 cM against TeoNAM's
#                         intended ~0.2 cM, because 64 junctions saturated at exactly
#                         1001.506 cM (rf = 0.5) dominate the upper tail and distort the
#                         ecdf/spline fit. At 75 cM the singleton pass flags nothing.
# Fixed values instead of the data-driven threshold, because the gap distribution is
# still polluted by saturation and cannot calibrate itself.
FINE_THR   <- 10    # singleton isolation, cM -- above our median (2.3), near q95 (14.4)
ISLAND_GAP_CM <- 20 # cluster isolation, cM -- set by FVRZ 2026-07-31 from the sweep
                    # (agent/sweep_quirky.log). 10 left a 55.4 cM hole on chr5 where
                    # find_quirky had removed 3 markers that were filling it; 20 is the
                    # only setting that clears the >25 cM gaps without collapsing
                    # coverage. fine_thr is inert at >= 10 -- island_thr is the live one.
ISLAND_MAX_N  <- 5L # max markers in a flagged cluster (was 20)
LEN_LO <- 1348; LEN_HI <- 1596
# Markers exempted from the distortion filter, the LD prune and find_quirky, so a
# specific set can be forced into the map to see what it does. Empty by default -- set
# via the environment, e.g.
#   PROTECT_MARKERS="m1,m2" Rscript scripts/qtl/build_map.R
# This is a DIAGNOSTIC switch, not part of the pipeline: anything protected bypasses
# every QC test, so a map built with it set is not a QC'd map.
PROTECT <- setdiff(strsplit(Sys.getenv("PROTECT_MARKERS", ""), ",")[[1]], "")
NCLUST <- min(8L, parallel::detectCores())
dir.create("output/qtl", showWarnings = FALSE, recursive = TRUE)
source(MAPTOOLS); stopifnot(exists("find_quirky"))

# verbatim from teonam_qtl_permap.R:55-74
renorm_z <- function(x) {
  z <- rep(NA_real_, length(x)); ok <- is.finite(x)
  if (sum(ok) < 5L) return(z)
  d <- ecdf(x[ok])
  u <- suppressWarnings(predict(smooth.spline(x[ok], d(x[ok])), x[ok])$y)
  u <- pmin(pmax(u, 1e-6), 1 - 1 / sum(ok))
  z[ok] <- qnorm(u); z
}
is_outlier <- function(x) { z <- renorm_z(x); !is.na(z) & z > 1.96 }
chr_len <- function(m) sapply(m, function(v) max(v) - min(v))

rule("0. INPUT")
cr0 <- read.cross(format = "csvsr", dir = "",
                  genfile = file.path(D, "rqtl_gen_abh_all.csv"),
                  phefile = file.path(D, "rqtl_phe.csv"),
                  genotypes = c("A","H","B"), na.strings = c("-","NA"),
                  crosstype = "f2", estimate.map = FALSE)
n_in <- totmar(cr0)
cat("markers:", n_in, " individuals:", nind(cr0), "\n")
info <- fread(file.path(D, "abh_all_marker_info.tsv"),
              colClasses = list(character = "chr_v5"))

rule("0b. DROP INVERTED MARKERS (A/B label backwards)")
# WHY. est.rf/checkAlleles flags markers whose alleles appear switched. A marker whose
# A/B label is backwards is anticorrelated with BOTH physical neighbours: every
# individual reads as a double recombinant there, rf -> ~0.5 on both sides, and Haldane
# blows up. Flipping is exact -- cor(2 - X, Y) = -cor(X, Y) -- so a true inversion has
# the SAME |r| as a normal marker, only the sign differs.
#
# EVIDENCE that these are mislabelled rather than bad (agent/diagnose_switched_alleles.log):
#   median |r| with the weaker neighbour 0.82 vs a 0.92 genome-wide baseline -- as tightly
#   linked as normal markers; call rate 185/186 vs 186; MAF 0.481 vs 0.478.
#   All 27 have cml_class == "ref" (baseline 984 ref / 881 alt, p ~ 3e-8), and the cause is
#   NOT the allele call: cml_base == ref_refor for all 984 ref markers, tag_base ==
#   cml_base for all 984, no strand enrichment (inverted are LESS often reverse-strand),
#   no A/T-C/G ambiguity enrichment (p = 1). See agent/localise_ref_inversion.log.
#   The mechanism is unresolved; the correction is applied on the linkage evidence alone.
#
# Detection is done here rather than read from a file so the step is reproducible and the
# flipped markers are named in this log.
G0 <- pull.geno(cr0) - 1L                       # A/H/B codes 1/2/3 -> dosage 0/1/2
storage.mode(G0) <- "double"
mi0 <- data.table(marker = colnames(G0),
                  chr = info$cml_chr[match(colnames(G0), info$marker)],
                  pos = info$ref_pos[match(colnames(G0), info$marker)])
stopifnot(!any(is.na(mi0$chr)), !any(is.na(mi0$pos)))
mi0 <- mi0[order(chr, pos)]
G0 <- G0[, mi0$marker, drop = FALSE]
pair_r <- function(a, b) {
  ok <- !is.na(a) & !is.na(b)
  if (sum(ok) < 4L || sd(a[ok]) == 0 || sd(b[ok]) == 0) return(NA_real_)
  cor(a[ok], b[ok])
}
mi0[, `:=`(r_left = NA_real_, r_right = NA_real_)]
for (ch in unique(mi0$chr)) {
  idx <- which(mi0$chr == ch)
  for (k in seq_along(idx)) {
    i <- idx[k]
    if (k > 1L)          set(mi0, i, "r_left",  pair_r(G0[, i], G0[, idx[k - 1L]]))
    if (k < length(idx)) set(mi0, i, "r_right", pair_r(G0[, i], G0[, idx[k + 1L]]))
  }
}
inverted <- mi0[r_left < 0 & r_right < 0, marker]
cat(sprintf("markers anticorrelated with BOTH neighbours: %d of %d (%.2f%%)\n",
            length(inverted), nrow(mi0), 100 * length(inverted) / nrow(mi0)))
print(mi0[marker %in% inverted, .(marker, chr, pos_Mb = round(pos / 1e6, 1),
                                  r_left = round(r_left, 3), r_right = round(r_right, 3))])
print(info[marker %in% inverted, .N, by = cml_class])
# REMOVED, not flipped -- FVRZ's call, 2026-07-31. Flipping would rescue the marker on
# the assumption that only its label is wrong; dropping makes no such assumption. The
# alignment filter at encoding already removed 22 of the original 27 on independent
# evidence (MAPQ/XA), so what remains is a small residue and the cost of dropping it is
# low. Note this also breaks A = CML530 in the other direction: a flipped marker would
# have kept its genotypes with A meaning the tester allele, which is why it is not done.
if (length(inverted)) {
  cr0 <- drop.markers(cr0, inverted)
  cat(sprintf("dropped %d inverted markers; %d remain\n", length(inverted), totmar(cr0)))
} else cat("no inverted markers to drop\n")
n_inverted <- length(inverted)

rule("1. SEGREGATION DISTORTION (per-chromosome renorm_z outliers, z > 1.96)")
G <- pull.geno(cr0); cv <- rep(names(cr0$geno), nmar(cr0)); nm <- colnames(G)
cnt <- cbind(colSums(G == 1, na.rm = TRUE),
             colSums(G == 2, na.rm = TRUE),
             colSums(G == 3, na.rm = TRUE))
tot <- rowSums(cnt); e <- tot %o% c(0.25, 0.50, 0.25)
chi <- rowSums((cnt - e)^2 / e)
cat("chi2 vs 1:2:1 -- median", round(median(chi),2),
    " q95", round(quantile(chi,.95),2), " max", round(max(chi),1), "\n")
dist_out <- logical(length(nm))
for (ch in unique(cv)) { i <- which(cv == ch); dist_out[i] <- is_outlier(chi[i]) }
cat("flagged:", sum(dist_out), sprintf("(%.1f%%)\n", 100*mean(dist_out)))
print(table(cv[dist_out]))
if (length(PROTECT)) dist_out[nm %in% PROTECT] <- FALSE   # diagnostic exemption
keep1 <- nm[!dist_out]
cr1 <- pull.markers(cr0, keep1)
cat("retained:", totmar(cr1), "\n")

rule(sprintf("1b. LD PRUNE (Pearson r >= %.2f between retained neighbours)", PRUNE_R))
# WHY. With 186 F2s there are 2 x 186 = 372 meioses, so the smallest estimable
# recombination fraction is 1/372 -> 0.27 cM. A single miscalled genotype inside a
# haplotype block reads as a DOUBLE crossover, so each error adds 2/372 = 0.54 cM of
# spurious length. Spurious length therefore scales with MARKER COUNT, not with
# error.prob -- which is why the error.prob sweep never plateaued. Pruning removes
# genotypes that contribute error without contributing resolution.
#
# findDupMarkers was tried here first and is NOT sufficient: it requires IDENTICAL
# genotype patterns and found only 21 groups / 28 markers / 1.5% of genotypes. Markers
# that differ at a few individuals are not caught, so the criterion has to be
# correlation, not identity.
#
# THE CRITERION. Teng & Xu 2026 (TAG, PMC12906585) Eq. 9: in an F_t (t >= 2) from inbred
# founders, with dosage coding X in {0,1,2}, the Pearson correlation between two loci is
#     r = 1 - 2*theta
# so with Haldane (their Eq. 29), since 1 - 2*theta = r,
#     d_cM = -50 * ln(r)
# This gives a genetic axis from correlation alone -- no HMM -- which is what makes
# pruning BEFORE est.map possible. A = CML530 on every marker with no phase flips, so
# the coding is coupling-consistent and the SIGN of r is meaningful. Negative r is set
# to 0 (Teng & Xu, Discussion), giving theta = 0.5.
#
# THE ALGORITHM. Keeping markers pairwise >= d apart is Maximum Distance-d Independent
# Set. On a general graph it is NP-hard (Jena, Jallu, Das & Nandy 2018, FAW 2018), which
# is why a solver over an O(n^2) matrix is heavy and only heuristic. On a 1-D axis the
# markers form an interval family and it collapses to the interval-scheduling greedy --
# sort, keep the leftmost, skip everything within d, repeat -- which is O(n) and EXACT.
# Taken from zealhmm/scripts/teonam_gwas118k_thin01.R:38-51.
#
# Runs AFTER the distortion filter, per the pipeline order, so a distorted marker cannot
# be retained as the representative of a region.
Gd <- pull.geno(cr1) - 1L                 # A/H/B codes 1/2/3 -> dosage 0/1/2
storage.mode(Gd) <- "double"
mi <- info[match(colnames(Gd), marker), .(marker, chr = cml_chr, pos = ref_pos)]
mi <- mi[order(chr, pos)]                 # fixed marker order from CML530
Gd <- Gd[, mi$marker, drop = FALSE]

adj_r <- function(cols) vapply(seq_len(length(cols) - 1L), function(k) {
  a <- Gd[, cols[k]]; b <- Gd[, cols[k + 1L]]
  ok <- !is.na(a) & !is.na(b)
  if (sum(ok) < 4L) return(NA_real_)
  if (sd(a[ok]) == 0 || sd(b[ok]) == 0) return(NA_real_)
  cor(a[ok], b[ok])
}, numeric(1))

axis_cm <- function(m) {
  r  <- adj_r(m$marker)
  th <- pmin((1 - pmax(r, 0)) / 2, 0.49)  # cap keeps Haldane finite
  th[is.na(th)] <- 0.49
  c(0, cumsum(-50 * log(1 - 2 * th)))
}
keep_thin <- unlist(lapply(split(mi, mi$chr), function(m) {
  cm <- axis_cm(m)
  last <- -Inf; k <- logical(length(cm))
  for (i in seq_along(cm)) if (cm[i] - last >= THIN_CM) { k[i] <- TRUE; last <- cm[i] }
  m$marker[k]
}), use.names = FALSE)

if (length(PROTECT)) keep_thin <- union(keep_thin, intersect(PROTECT, colnames(Gd)))
n_pruned <- totmar(cr1) - length(keep_thin)
cat(sprintf("r = %.3f  ->  minimum spacing %.4f cM on the correlation axis\n",
            PRUNE_R, THIN_CM))
cat(sprintf("markers %d -> %d (pruned %d, %.1f%%)\n", totmar(cr1), length(keep_thin),
            n_pruned, 100 * n_pruned / totmar(cr1)))
cr1 <- pull.markers(cr1, keep_thin)
cat("per chromosome:\n"); print(nmar(cr1))

rule(sprintf("2. ROUND 1 est.map (PRELIMINARY, error.prob=%g, %s)", ERROR_PROB, MAP_FUN))
t0 <- Sys.time()
m1 <- est.map(cr1, error.prob = ERROR_PROB, map.function = MAP_FUN,
              maxit = 10000, tol = 1e-6, n.cluster = NCLUST)
cat(sprintf("elapsed %.1f min\n", as.numeric(difftime(Sys.time(), t0, units="mins"))))
l1 <- chr_len(m1); print(round(l1, 1))
cat(sprintf("round-1 TOTAL: %.1f cM over %d markers\n", sum(l1), totmar(cr1)))
cat(sprintf("expected spacing if the map were 1800 cM: %.3f cM\n", 1800/totmar(cr1)))
g1 <- unlist(lapply(m1, function(v) diff(v[order(v)])), use.names = FALSE)
cat(sprintf("gaps: median %.3f  q95 %.3f  q99 %.3f  max %.3f cM | >10cM: %d  >25cM: %d\n",
            median(g1), quantile(g1,.95), quantile(g1,.99), max(g1),
            sum(g1 > 10), sum(g1 > 25)))

rule("3. find_quirky (data-driven fine threshold + island rule)")
gap_out_thr <- { o <- is_outlier(g1); if (any(o)) min(g1[o]) else Inf }
cat(sprintf("data-driven threshold would be %.3f cM -- REPORTED ONLY, not used\n",
            gap_out_thr))
cat(sprintf("using fine_thr = %g, island_thr = %g, island_max_n = %d\n",
            FINE_THR, ISLAND_GAP_CM, ISLAND_MAX_N))
quirky <- unlist(lapply(m1, find_quirky, fine_thr = FINE_THR,
                        island_thr = ISLAND_GAP_CM, island_max_n = ISLAND_MAX_N),
                 use.names = FALSE)
if (length(PROTECT)) quirky <- setdiff(quirky, PROTECT)   # diagnostic exemption
cat("flagged:", length(quirky), sprintf("(%.1f%% of round-1 set)\n",
    100*length(quirky)/totmar(cr1)))
if (length(quirky)) print(table(info[marker %in% quirky, cml_chr]))
# keep_thin, NOT keep1: keep1 still contains the pruned markers, so building keep2 from
# it would reinstate every pruned marker for round 2 and cancel the prune entirely.
keep2 <- setdiff(keep_thin, quirky)
cr2 <- pull.markers(cr0, keep2)
cat("retained:", totmar(cr2), "\n")

rule("4. ROUND 2 est.map (REFINED)")
t0 <- Sys.time()
m2 <- est.map(cr2, error.prob = ERROR_PROB, map.function = MAP_FUN,
              maxit = 10000, tol = 1e-6, n.cluster = NCLUST)
cat(sprintf("elapsed %.1f min\n", as.numeric(difftime(Sys.time(), t0, units="mins"))))
cr2 <- replace.map(cr2, m2)
l2 <- chr_len(m2)
# PHYSICAL COVERAGE, not just cM. A chromosome can look monotonic and sigmoid while
# retaining only a fragment of its length -- earlier builds lost whole arms (chr7 kept
# 11% of its length, chr10 15%) and I reported only total cM, which conceals that.
# Coordinates are CML530's now, so coverage must be against CML530 chromosome lengths,
# not B73's. Read from the BAM header so they cannot drift from the assembly the
# positions were taken from.
CML_LEN <- local({
  sq <- system(sprintf("samtools view -H %s | awk '$1==\"@SQ\"{print $2\"\\t\"$3}'",
                       file.path(D, "CML530", "tags_vs_CML530.bam")), intern = TRUE)
  nm  <- sub("^SN:", "", sub("\t.*$", "", sq))
  len <- as.numeric(sub("^LN:", "", sub("^.*\t", "", sq)))
  names(len) <- sub("^chr0?", "", nm)
  len[as.character(1:10)]
})
stopifnot(!any(is.na(CML_LEN)))
inmap <- lapply(m2, names)
pos_of <- function(ch) info[marker %in% inmap[[ch]], sort(ref_pos)]
tab <- data.table(chr = names(l2), markers = nmar(cr2), cM = round(l2, 1))
tab[, mean_gap_cM := round(cM / (markers - 1), 3)]
tab[, first_Mb := round(sapply(chr, function(ch) min(pos_of(ch))) / 1e6, 1)]
tab[, last_Mb  := round(sapply(chr, function(ch) max(pos_of(ch))) / 1e6, 1)]
tab[, chr_Mb   := round(CML_LEN[chr] / 1e6, 1)]
tab[, span_Mb  := round(last_Mb - first_Mb, 1)]
tab[, pct_covered := round(100 * span_Mb / chr_Mb, 1)]
tab[, lost_start_Mb := first_Mb]
tab[, lost_end_Mb   := round(chr_Mb - last_Mb, 1)]
setcolorder(tab, c("chr","markers","cM","mean_gap_cM","first_Mb","last_Mb",
                   "chr_Mb","pct_covered","lost_start_Mb","lost_end_Mb"))
print(tab[order(as.integer(chr))])
cat(sprintf("\ngenome-wide physical coverage: %.1f%% (%.0f of %.0f Mb)\n",
            100*sum(tab$span_Mb)/sum(tab$chr_Mb), sum(tab$span_Mb), sum(tab$chr_Mb)))
cat("lost_start_Mb / lost_end_Mb are UNMAPPED regions -- a QTL there cannot be found.\n")
cat(sprintf("\nTOTAL: %.1f cM over %d markers\n", sum(l2), totmar(cr2)))
g2 <- unlist(lapply(m2, function(v) diff(v[order(v)])), use.names = FALSE)
cat(sprintf("gaps: median %.3f  q95 %.3f  q99 %.3f  max %.3f cM | >10cM: %d  >25cM: %d\n",
            median(g2), quantile(g2,.95), quantile(g2,.99), max(g2),
            sum(g2 > 10), sum(g2 > 25)))

rule("5. GUARDRAILS")
dr <- (n_in - totmar(cr2)) / n_in
cat(sprintf("markers %d -> %d (inverted -%d, distortion -%d, LD prune -%d, quirky -%d, drop rate %.1f%%)\n",
            n_in, totmar(cr2), n_inverted, sum(dist_out), n_pruned, length(quirky), 100*dr))
if (dr > 0.10) cat("WARN: drop rate exceeds 10%\n")
cat(sprintf("length %.1f cM vs %d-%d cM: %s\n", sum(l2), LEN_LO, LEN_HI,
            ifelse(sum(l2) >= LEN_LO & sum(l2) <= LEN_HI, "PASS", "OUT OF RANGE")))

rule("6. OUTPUT + PLOT")
saveRDS(cr2, file.path(D, "rqtl_cross_map_teonamqc.rds"))
qc <- info[, .(marker, cml_chr, ref_pos, chr_v5, pos_v5, cml_class, mapq, maf, p121)]
qc[, `:=`(dropped_distortion = marker %in% nm[dist_out],
          dropped_quirky = marker %in% quirky,
          in_map = marker %in% unlist(lapply(m2, names)))]
fwrite(qc, file.path(D, "map_qc_teonamqc.tsv"), sep = "\t")

# PUBLIC MAP. data/qtl/ is gitignored (raw DArTseq genotypes are not redistributable),
# so the map itself -- marker names, linkage group, cM, and both coordinate systems --
# is written to a TRACKED path. It contains no individual genotypes.
pub <- rbindlist(lapply(names(m2), function(ch) {
  v <- sort(m2[[ch]])
  data.table(marker = names(v), lg = as.integer(ch),
             cM = round(as.numeric(v) - min(as.numeric(v)), 4))
}))
pub <- merge(pub, info[, .(marker, cml_chr, cml_pos = ref_pos,
                           b73v5_chr = chr_v5, b73v5_pos = pos_v5,
                           cml_allele = cml_class, maf, p121)],
             by = "marker", sort = FALSE)
setorder(pub, lg, cM)
setcolorder(pub, c("marker", "lg", "cM", "cml_chr", "cml_pos",
                   "b73v5_chr", "b73v5_pos", "cml_allele", "maf", "p121"))
fwrite(pub, "data/f2_genetic_map.tsv", sep = "\t")
cat(sprintf("wrote data/f2_genetic_map.tsv  (%d markers, public)\n", nrow(pub)))

mkp <- rbindlist(lapply(names(m2), function(ch)
  data.table(chr = ch, cm = m2[[ch]] - min(m2[[ch]]))))
mkp[, chr := factor(chr, levels = as.character(1:10))]
p <- ggplot(mkp, aes(x = chr, y = cm)) +
  # Markers are drawn OPAQUE BLACK as horizontal ticks. The TeoNAM version used
  # grey30 at alpha 0.20, which reads as a density on a 9,000-marker map but makes
  # individual markers invisible on a 1,152-marker one.
  geom_segment(aes(xend = chr, y = 0, yend = cm), linewidth = 9, colour = "grey88") +
  geom_segment(aes(x = as.numeric(chr) - 0.26, xend = as.numeric(chr) + 0.26,
                   y = cm, yend = cm), linewidth = 0.35, colour = "black") +
  scale_y_reverse(expand = expansion(mult = c(0.02, 0.04))) +
  labs(x = "chromosome (CML530)", y = "position (cM)",
       title = sprintf("F2 map -- FVRZ QC on %d encoded markers: %d markers, %.1f cM", n_in,
                       totmar(cr2), sum(l2))) +
  theme_classic(base_size = 13)
ggsave("output/qtl/f2_map_teonamqc_ideogram.png", p, width = 9, height = 6, dpi = 150, bg = "white")

mk2 <- rbindlist(lapply(names(m2), function(ch)
  data.table(chr = ch, marker = names(m2[[ch]]),
             cm = as.numeric(m2[[ch]]) - min(m2[[ch]]))))
mk2 <- merge(mk2, info[, .(marker, ref_pos)], by = "marker")
mk2[, `:=`(mb = ref_pos/1e6, chr = factor(chr, levels = as.character(1:10)))]
q <- ggplot(mk2[order(chr, mb)], aes(mb, cm)) +
  geom_line(colour = "grey60", linewidth = 0.4) +
  geom_point(size = 0.4, colour = "grey20") +
  facet_wrap(~chr, scales = "free", ncol = 5) +
  labs(x = "physical position (Mb, CML530)", y = "genetic position (cM)",
       title = "Marey map -- FVRZ QC") +
  theme_classic(base_size = 11)
ggsave("output/qtl/f2_map_teonamqc_marey.png", q, width = 12, height = 5.5, dpi = 150, bg = "white")
cat("wrote output/qtl/f2_map_teonamqc_{ideogram,marey}.png\n")
