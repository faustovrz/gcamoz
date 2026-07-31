#!/usr/bin/env Rscript
# qtl_build_map.R -- estimate the F2 genetic map on the FIXED B73 v5 physical order.
#
# APPROACH (agreed): do NOT order markers de novo. With n=166 the ordering problem
# is under-determined and would produce a worse map than the physical order. Fix
# the order from B73 v5, estimate DISTANCES with est.map, and use est.rf only to
# FLAG markers whose recombination pattern contradicts that order.
#
# PIPELINE -- the two-round est.map bootstrap from zealhmm/scripts/teonam_qtl_permap.R,
# adapted from BC1S4 to F2. Bad markers cannot be found on the map you are trying
# to build, so build a throwaway map first, use it to find them, then rebuild:
#
#   distortion drop -> ROUND 1 est.map (preliminary) -> quirky drop -> ROUND 2 (refined)
#
# ADAPTATIONS FROM THE TeoNAM VERSION
#   * expected genotype frequencies: BC1S4 (0.734, 0.031, 0.234) -> F2 (0.25, 0.50, 0.25)
#   * map.function: haldane -> kosambi (maize convention)
#   * error.prob: 0.001 -> 0.01. DArT read depth is modest (median 6.9x ref /
#     5.5x alt); the default 1e-4 would force real genotyping error into map length.
#   * length guardrail: Chen's 1348-1596 cM (TeoNAM subpop) -> ~1400-1800 cM (maize F2)
#
# find_quirky()/find_quirky_islands() are SOURCED from zealhmm, not reimplemented,
# so the island rule stays identical across both projects.
#
# A NOTE ON THE DISTORTION FILTER: is_outlier() flags a fixed ~2.5% upper tail by
# construction (z > 1.96), so it ALWAYS removes roughly that fraction whether or
# not anything is wrong. A low drop rate is therefore NOT validation. Total map
# length against the maize expectation is what validates the build.
#
# Reads   data/qtl/derived/rqtl_gen_abh.csv, rqtl_phe.csv, abh_marker_info.tsv
# Writes  data/qtl/derived/rqtl_cross_map.rds   (cross object with the refined map)
#         data/qtl/derived/map_qc.tsv           (per-marker drop reasons)
# Usage: Rscript scripts/qtl/qtl_build_map.R > agent/qtl_build_map.log 2>&1

suppressMessages({library(qtl); library(data.table)})
rule <- function(x) cat("\n", strrep("=", 70), "\n", x, "\n", strrep("=", 70), "\n", sep = "")

D        <- "data/qtl/derived"
# Prefer the LINKAGE-PHASED genotypes. Unphased markers make rf -> 1-rf and blew
# the first build up to 161,170 cM; qtl_phase_markers.R took phase consistency
# from 92.4% to 99.3%.
GENFILE  <- if (file.exists(file.path(D, "rqtl_gen_abh_phased.csv")))
              file.path(D, "rqtl_gen_abh_phased.csv") else
              file.path(D, "rqtl_gen_abh.csv")
PHEFILE  <- file.path(D, "rqtl_phe.csv")
MKINFO   <- file.path(D, "abh_marker_info.tsv")
OUT_RDS  <- file.path(D, "rqtl_cross_map.rds")
OUT_QC   <- file.path(D, "map_qc.tsv")
MAPTOOLS <- "/Users/fvrodriguez/repos/zealhmm/scripts/map_tools.R"

ERROR_PROB    <- 0.01
MAP_FUN       <- "kosambi"
ISLAND_GAP_CM <- 2
ISLAND_MAX_N  <- 20L
LEN_LO        <- 1400
LEN_HI        <- 1800
NCLUST        <- min(8L, parallel::detectCores())

rule("0. PREFLIGHT")
stopifnot(file.exists(GENFILE), file.exists(PHEFILE), file.exists(MKINFO))
if (!file.exists(MAPTOOLS)) stop("map_tools.R not found at ", MAPTOOLS)
source(MAPTOOLS)   # find_quirky(), find_quirky_islands()
stopifnot(exists("find_quirky"), exists("find_quirky_islands"))
cat("sourced find_quirky from", MAPTOOLS, "\n")
cat(sprintf("error.prob=%.3f  map.function=%s  island_thr=%g cM  island_max_n=%d\n",
            ERROR_PROB, MAP_FUN, ISLAND_GAP_CM, ISLAND_MAX_N))

# --- airmine empirical-CDF -> qnorm renormalisation; upper-tail outliers z>1.96 -
# Vendored from zealhmm/scripts/teonam_qtl_permap.R (lines 55-74) rather than
# sourced, because that file is a run script, not a function library.
renorm_z <- function(x) {
  z <- rep(NA_real_, length(x)); ok <- is.finite(x)
  if (sum(ok) < 5L) return(z)
  d <- ecdf(x[ok])
  u <- suppressWarnings(predict(smooth.spline(x[ok], d(x[ok])), x[ok])$y)
  u <- pmin(pmax(u, 1e-6), 1 - 1 / sum(ok))
  z[ok] <- qnorm(u); z
}
is_outlier <- function(x) { z <- renorm_z(x); !is.na(z) & z > 1.96 }

rule("1. READ CROSS")
cr0 <- read.cross(format = "csvsr", dir = "", genfile = GENFILE, phefile = PHEFILE,
                  genotypes = c("A", "H", "B"), na.strings = c("-", "NA"),
                  crosstype = "f2", estimate.map = FALSE)
cat("individuals:", nind(cr0), " markers:", totmar(cr0), " chromosomes:", nchr(cr0), "\n")
mkinfo <- fread(MKINFO, colClasses = list(character = "chr_v5"))
cat("marker info rows:", nrow(mkinfo), "\n")
cat("polarity_source:\n"); print(mkinfo[, .N, by = polarity_source])

rule("2. SEGREGATION DISTORTION (F2 expectation 0.25 / 0.50 / 0.25)")
G <- pull.geno(cr0)                       # individuals x markers, 1=AA 2=AB 3=BB
mk_chr <- rep(names(cr0$geno), nmar(cr0))
mk_nm  <- colnames(G)
cnt <- cbind(AA = colSums(G == 1, na.rm = TRUE),
             AB = colSums(G == 2, na.rm = TRUE),
             BB = colSums(G == 3, na.rm = TRUE))
tot <- rowSums(cnt)
f2_exp <- c(0.25, 0.50, 0.25)
e <- tot %o% f2_exp
chi <- rowSums((cnt - e)^2 / e)
# per-chromosome outlier detection: absorbs genuine chromosome-wide distortion
# (a segregating lethal, an inversion) instead of deleting the chromosome.
dist_out <- logical(length(mk_nm))
for (ch in unique(mk_chr)) {
  i <- which(mk_chr == ch)
  dist_out[i] <- is_outlier(chi[i])
}
cat(sprintf("chi2 vs 1:2:1 -- median %.2f  q95 %.2f  max %.2f\n",
            median(chi), quantile(chi, .95), max(chi)))
cat("flagged as distorted:", sum(dist_out), sprintf("(%.1f%%)\n",
    100 * mean(dist_out)))
cat("per chromosome:\n"); print(table(mk_chr[dist_out]))
cat("\nNB: this rule removes a ~2.5% upper tail BY CONSTRUCTION. A low drop rate\n")
cat("is not evidence the data are clean; map length is the real check.\n")
keep1 <- mk_nm[!dist_out]
cr1 <- pull.markers(cr0, keep1)
cat("markers after distortion drop:", totmar(cr1), "\n")

rule("2b. LOCAL-LINKAGE FILTER (drop markers unlinked to their neighbours)")
# WHY THIS IS NEEDED AND WHY IT COMES FIRST.
# Phase correction fixed the SIGN of adjacent correlations (92.4% -> 99.3%), but
# sign is not strength: ~5% of adjacent pairs have |cor| < 0.3, i.e. no linkage at
# ANY polarity. Those are paralogous or erroneous markers, and each one injects
# ~50 cM of fictional distance. Left in, they inflated round 1 to 12,865 cM, which
# in turn drove find_quirky's data-driven gap threshold to an absurd 91.8 cM and
# made it flag 623 markers (45%). Garbage in, garbage out.
#
# So remove them DIRECTLY on an interpretable criterion -- a marker must correlate
# with at least one physical neighbour -- rather than letting est.map produce a
# nonsense map and asking the island rule to clean up after it.
LINK_MIN <- 0.30      # |cor| a marker must reach with at least one neighbour
LINK_K   <- 5L        # neighbours each side to consider
G1 <- pull.geno(cr1); c1 <- rep(names(cr1$geno), nmar(cr1)); d1 <- G1 - 1L
best <- rep(NA_real_, ncol(d1))
for (ch in unique(c1)) {
  idx <- which(c1 == ch); n <- length(idx)
  if (n < 2) next
  dd <- d1[, idx, drop = FALSE]
  for (k in seq_len(n)) {
    nb <- setdiff(max(1, k - LINK_K):min(n, k + LINK_K), k)
    cs <- sapply(nb, function(j) {
      a <- dd[, k]; b <- dd[, j]; ok <- is.finite(a) & is.finite(b)
      if (sum(ok) < 30) return(NA_real_)
      abs(cor(a[ok], b[ok]))
    })
    best[idx[k]] <- suppressWarnings(max(cs, na.rm = TRUE))
  }
}
cat("max |cor| with any neighbour (per marker):\n")
print(round(summary(best), 3))
unlinked <- colnames(d1)[is.finite(best) & best < LINK_MIN]
cat(sprintf("\nmarkers below |cor| %.2f with EVERY neighbour: %d (%.1f%%)\n",
            LINK_MIN, length(unlinked), 100 * length(unlinked) / ncol(d1)))
cat("per chromosome:\n"); print(table(c1[colnames(d1) %in% unlinked]))
keep1b <- setdiff(keep1, unlinked)
cr1 <- pull.markers(cr0, keep1b)
cat("markers after distortion + local-linkage drop:", totmar(cr1), "\n")

rule("3. ROUND 1 -- PRELIMINARY est.map (locates quirky markers only)")
t0 <- Sys.time()
m1 <- est.map(cr1, error.prob = ERROR_PROB, map.function = MAP_FUN,
              maxit = 10000, tol = 1e-6, n.cluster = NCLUST)
cat(sprintf("elapsed %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
len1 <- sapply(m1, function(v) max(v) - min(v))
cat("round-1 length per chromosome (cM):\n"); print(round(len1, 1))
cat(sprintf("round-1 TOTAL: %.1f cM\n", sum(len1)))

rule("4. QUIRKY MARKERS (data-driven gap threshold + isolated-cluster rule)")
gaps_all <- unlist(lapply(m1, function(v) diff(v[order(v)])), use.names = FALSE)
gap_out_thr <- { o <- is_outlier(gaps_all); if (any(o)) min(gaps_all[o]) else Inf }
cat(sprintf("inter-marker gaps: median %.3f  q99 %.3f  max %.3f cM\n",
            median(gaps_all), quantile(gaps_all, .99), max(gaps_all)))
cat(sprintf("data-driven fine gap threshold (smallest OUTLIER gap): %.3f cM\n", gap_out_thr))
quirky <- unlist(lapply(m1, find_quirky, fine_thr = gap_out_thr,
                        island_thr = ISLAND_GAP_CM, island_max_n = ISLAND_MAX_N),
                 use.names = FALSE)
cat("quirky markers flagged:", length(quirky), "\n")
if (length(quirky)) {
  qi <- mkinfo[marker %in% quirky]
  cat("per chromosome:\n"); print(qi[, .N, by = chr_v5][order(as.integer(chr_v5))])
}
# keep1b, NOT keep1: keep1 predates the local-linkage filter, so using it here
# would silently re-add every unlinked marker that section 2b removed.
keep2 <- setdiff(keep1b, quirky)
cr2 <- pull.markers(cr0, keep2)
cat("markers after distortion + quirky drop:", totmar(cr2), "\n")

rule("5. ROUND 2 -- REFINED est.map")
t0 <- Sys.time()
m2 <- est.map(cr2, error.prob = ERROR_PROB, map.function = MAP_FUN,
              maxit = 10000, tol = 1e-6, n.cluster = NCLUST)
cat(sprintf("elapsed %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
cr2 <- replace.map(cr2, m2)
len2 <- sapply(m2, function(v) max(v) - min(v))
tab <- data.table(chr = names(len2), n_markers = nmar(cr2),
                  cM = round(len2, 1),
                  Mb = round(sapply(names(len2), function(ch)
                    diff(range(mkinfo[marker %in% colnames(pull.geno(cr2, chr = ch)),
                                      pos_v5])) / 1e6), 1))
tab[, cM_per_Mb := round(cM / Mb, 2)]
print(tab)
cat(sprintf("\nREFINED TOTAL: %.1f cM over %d markers\n", sum(len2), totmar(cr2)))

rule("6. GUARDRAILS")
n_in <- totmar(cr0); drop_rate <- (n_in - totmar(cr2)) / n_in
cat(sprintf("markers %d -> %d  (distortion -%d, quirky -%d, drop rate %.1f%%)\n",
            n_in, totmar(cr2), sum(dist_out), length(quirky), 100 * drop_rate))
if (drop_rate > 0.10) cat("WARN: drop rate exceeds 10%\n")
cat(sprintf("total length %.1f cM against the maize F2 expectation %d-%d cM: %s\n",
            sum(len2), LEN_LO, LEN_HI,
            ifelse(sum(len2) >= LEN_LO & sum(len2) <= LEN_HI, "PASS",
                   "OUT OF RANGE -- see note")))
if (sum(len2) > LEN_HI)
  cat("  A length far above expectation means residual genotyping error is\n",
      "  inflating distances: raise error.prob and re-run.\n")
if (sum(len2) < LEN_LO)
  cat("  A length far below expectation suggests over-aggressive filtering or\n",
      "  an error.prob absorbing real recombination.\n")

rule("7. est.rf -- does recombination contradict the physical order?")
cr2 <- est.rf(cr2)
rf <- pull.rf(cr2, what = "lod")
worst <- rbindlist(lapply(names(cr2$geno), function(ch) {
  mk <- colnames(pull.geno(cr2, chr = ch)); if (length(mk) < 3) return(NULL)
  sub <- rf[mk, mk]
  # LOD to the immediate physical neighbour: low = order likely wrong locally
  adj <- sapply(seq_len(length(mk) - 1), function(i) sub[i, i + 1])
  data.table(chr = ch, marker = mk[-length(mk)], lod_next = adj)
}))
cat("LOD to next physical marker -- low values flag local order conflict:\n")
print(summary(worst$lod_next))
cat("markers with lod_next < 3:", sum(worst$lod_next < 3, na.rm = TRUE), "\n")
cat("(NOT dropped here -- reported for inspection; the physical order stands)\n")

rule("8. OUTPUTS")
qc <- mkinfo[, .(marker, chr_v5, pos_v5, polarity_source, cml_class)]
qc[, `:=`(dropped_distortion = marker %in% mk_nm[dist_out],
          dropped_quirky     = marker %in% quirky,
          in_map             = marker %in% unlist(lapply(m2, names)))]
cm <- unlist(lapply(m2, function(v) v - min(v)))
qc[, cM := cm[match(marker, unlist(lapply(m2, names)))]]
qc <- merge(qc, worst[, .(marker, lod_next)], by = "marker", all.x = TRUE)
setorder(qc, chr_v5, pos_v5)
fwrite(qc, OUT_QC, sep = "\t")
saveRDS(cr2, OUT_RDS)
cat("wrote", OUT_QC, "rows:", nrow(qc), "\n")
cat("wrote", OUT_RDS, "\n")
print(summary(cr2))
