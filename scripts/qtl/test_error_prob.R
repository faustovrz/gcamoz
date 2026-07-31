#!/usr/bin/env Rscript
# test_error_prob.R -- sweep est.map's error.prob on the CLEAN marker set.
#
# WHY NOW AND NOT BEFORE. Every earlier attempt at this would have been fitting noise,
# because the marker set was contaminated: non-segregating markers, co-tag duplicates,
# markers contradicting the A = CML530 encoding, and -- worst -- find_quirky running with
# island_thr = 2 against a median gap of 2.3 cM, which deleted 79-91% of markers. That is
# fixed (fine_thr = 10, island_thr = 10, island_max_n = 5; drop rate now 8.7%, PASSES).
#
# THE REMAINING SYMPTOM IS UNIFORM, NOT LOCALISED:
#   total 4,736.6 cM vs 1,348-1,596 expected           ~3x too long
#   median gap 2.045 cM vs ~1.06 expected at 1,702 markers   ~2x too long
#   all ten chromosomes 242-794 cM, mean gaps 2.1-3.2 cM     evenly inflated
#   max gap 55.4 cM, only 3 gaps > 25 cM                     no saturation left
# A uniform inflation of every ordinary interval, with no chromosome or region carrying
# it, is what happens when est.map reads GENOTYPING ERROR as crossovers. error.prob is
# the parameter that controls exactly that.
#
# error.prob = 0.001 came from FVRZ's BC1S4 pipeline, where lines are near-homozygous
# and calls are easy. A 186-plant F2 on DArTseq at ~6x read depth is a different regime.
#
# HOW TO READ THIS -- AND THE TRAP.
# Do NOT pick the value that lands the total in 1,348-1,596; that is fitting the answer.
# Pick where the curve FLATTENS: length should be steeply sensitive while the model is
# absorbing real crossovers, and insensitive once it is only absorbing genuine error.
# If the curve never flattens inside a plausible range, error.prob is NOT the
# explanation and something else is.
#
# READ-ONLY with respect to the cross; writes a plot and a summary table.
# Usage: Rscript scripts/qtl/test_error_prob.R > agent/test_error_prob.log 2>&1

suppressMessages({library(qtl); library(data.table); library(ggplot2)})
rule <- function(x) cat("\n", strrep("=", 70), "\n", x, "\n", strrep("=", 70), "\n", sep = "")

D <- "data/qtl/derived"
MAPTOOLS <- "/Users/fvrodriguez/repos/zealhmm/scripts/map_tools.R"
MAP_FUN <- "haldane"
FINE_THR <- 10; ISLAND_GAP_CM <- 10; ISLAND_MAX_N <- 5L
EPS <- c(1e-4, 1e-3, 2.5e-3, 5e-3, 1e-2, 2e-2, 4e-2)
NCLUST <- min(8L, parallel::detectCores())
dir.create("output/qtl", showWarnings = FALSE, recursive = TRUE)
source(MAPTOOLS)

renorm_z <- function(x) {
  z <- rep(NA_real_, length(x)); ok <- is.finite(x)
  if (sum(ok) < 5L) return(z)
  d <- ecdf(x[ok]); u <- suppressWarnings(predict(smooth.spline(x[ok], d(x[ok])), x[ok])$y)
  z[ok] <- qnorm(pmin(pmax(u, 1e-6), 1 - 1/sum(ok))); z
}
is_outlier <- function(x) { z <- renorm_z(x); !is.na(z) & z > 1.96 }
chr_len <- function(m) sapply(m, function(v) max(v) - min(v))

rule("1. CROSS + FVRZ DISTORTION FILTER (identical to build_map.R)")
cr0 <- read.cross(format = "csvsr", dir = "",
                  genfile = file.path(D, "rqtl_gen_abh_all.csv"),
                  phefile = file.path(D, "rqtl_phe.csv"),
                  genotypes = c("A","H","B"), na.strings = c("-","NA"),
                  crosstype = "f2", estimate.map = FALSE)
G <- pull.geno(cr0); cv <- rep(names(cr0$geno), nmar(cr0)); nm <- colnames(G)
cnt <- cbind(colSums(G==1, na.rm=TRUE), colSums(G==2, na.rm=TRUE), colSums(G==3, na.rm=TRUE))
tot <- rowSums(cnt); e <- tot %o% c(.25,.5,.25); chi <- rowSums((cnt-e)^2/e)
dout <- logical(length(nm))
for (ch in unique(cv)) { i <- which(cv==ch); dout[i] <- is_outlier(chi[i]) }
cr1 <- pull.markers(cr0, nm[!dout])
cat("markers in:", totmar(cr0), " after distortion filter:", totmar(cr1), "\n")

rule("2. SWEEP -- full two-round pipeline at each error.prob")
res <- rbindlist(lapply(EPS, function(ep) {
  t0 <- Sys.time()
  m1 <- est.map(cr1, error.prob = ep, map.function = MAP_FUN,
                maxit = 10000, tol = 1e-6, n.cluster = NCLUST)
  g1 <- unlist(lapply(m1, function(v) diff(v[order(v)])), use.names = FALSE)
  qk <- unlist(lapply(m1, find_quirky, fine_thr = FINE_THR,
                      island_thr = ISLAND_GAP_CM, island_max_n = ISLAND_MAX_N),
               use.names = FALSE)
  cr2 <- pull.markers(cr0, setdiff(nm[!dout], qk))
  m2 <- est.map(cr2, error.prob = ep, map.function = MAP_FUN,
                maxit = 10000, tol = 1e-6, n.cluster = NCLUST)
  l2 <- chr_len(m2)
  g2 <- unlist(lapply(m2, function(v) diff(v[order(v)])), use.names = FALSE)
  cat(sprintf("error.prob %-7g r1 %8.1f  quirky %4d  markers %4d  r2 TOTAL %7.1f cM  median gap %.3f  (%.1f min)\n",
              ep, sum(chr_len(m1)), length(qk), totmar(cr2), sum(l2),
              median(g2), as.numeric(difftime(Sys.time(), t0, units="mins"))))
  data.table(error_prob = ep, round1_cM = sum(chr_len(m1)), n_quirky = length(qk),
             markers = totmar(cr2), total_cM = sum(l2),
             median_gap = median(g2), max_gap = max(g2),
             n_gap_gt10 = sum(g2 > 10), n_gap_gt25 = sum(g2 > 25))
}))

rule("3. SENSITIVITY -- WHERE DOES IT FLATTEN?")
res[, pct_change := c(NA, round(100*diff(total_cM)/head(total_cM,-1), 1))]
print(res)
cat("\nexpected median gap if the map were 1500 cM:",
    round(1500/median(res$markers), 3), "cM\n")
cat("target total: 1348-1596 cM (Chen). DO NOT pick to hit it -- pick the plateau.\n")

rule("4. PER-CHROMOSOME AT EACH error.prob")
det <- rbindlist(lapply(EPS, function(ep) {
  m1 <- est.map(cr1, error.prob = ep, map.function = MAP_FUN, maxit = 10000,
                tol = 1e-6, n.cluster = NCLUST)
  qk <- unlist(lapply(m1, find_quirky, fine_thr = FINE_THR,
                      island_thr = ISLAND_GAP_CM, island_max_n = ISLAND_MAX_N),
               use.names = FALSE)
  cr2 <- pull.markers(cr0, setdiff(nm[!dout], qk))
  m2 <- est.map(cr2, error.prob = ep, map.function = MAP_FUN, maxit = 10000,
                tol = 1e-6, n.cluster = NCLUST)
  data.table(error_prob = ep, chr = names(chr_len(m2)), cM = round(chr_len(m2), 1))
}))
print(dcast(det, chr ~ error_prob, value.var = "cM")[order(as.integer(chr))])

rule("5. PLOT")
p <- ggplot(res, aes(error_prob, total_cM)) +
  annotate("rect", xmin = min(EPS), xmax = max(EPS), ymin = 1348, ymax = 1596,
           fill = "grey85") +
  geom_line(colour = "grey50") + geom_point(size = 2) +
  scale_x_log10() +
  labs(x = "error.prob (log scale)", y = "total map length (cM)",
       title = "Map length vs assumed genotyping-error rate, clean marker set",
       subtitle = "grey band = Chen 1348-1596 cM; choose the plateau, not the band") +
  theme_classic(base_size = 12)
ggsave("output/qtl/f2_errorprob_sweep.png", p, width = 8, height = 5, dpi = 150, bg = "white")
fwrite(res, file.path(D, "errorprob_sweep.tsv"), sep = "\t")
cat("wrote output/qtl/f2_errorprob_sweep.png and", file.path(D, "errorprob_sweep.tsv"), "\n")
