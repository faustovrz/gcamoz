#!/usr/bin/env Rscript
# qtl_map_errorprob_sweep.R -- rebuild the map WITHOUT find_quirky and sweep error.prob.
#
# TWO CHANGES FROM qtl_build_map.R, both justified by the Marey plot
# (output/qtl/f2_map_marey.png), which is monotonic and sigmoid on all ten
# chromosomes -- i.e. marker ORDER and PHASE are correct:
#
# 1. DROP find_quirky ENTIRELY. Its purpose is catching markers whose position
#    contradicts the order. The Marey map shows none. It was dropping 407 markers
#    (33.8%) because it was written for the TeoNAM BC1S4 map, where the
#    99.99th-percentile gap is ~1 cM; here q99 is ~10 cM with wide pericentromeric
#    plateaus, so legitimate dense clumps separated by a few cM register as
#    "isolated islands" under island_thr = 2. I imported the rule without checking
#    that its assumptions transferred. They do not.
#
# 2. SWEEP error.prob instead of fixing it at 0.01. At 878 x 166 genotypes,
#    error.prob = 0.01 lets the HMM attribute ~1,450 calls to error, and at a
#    median gap of 0.39 cM each absorbed apparent double-crossover deletes two real
#    crossovers -- the classic cause of a too-short map (880 cM vs 1,400-1,800).
#
# HOW TO READ THE SWEEP -- AND THE TRAP TO AVOID.
# Do NOT pick the error.prob that makes the total match the literature; that is
# fitting the answer, which is the mistake I made earlier with the linkage
# threshold. Pick where the curve FLATTENS: length should be steeply sensitive to
# error.prob while the model is absorbing real crossovers, and insensitive once it
# is only absorbing genuine error. If the curve never flattens inside a plausible
# range, error.prob is not the explanation and something else is.
#
# Reads   data/qtl/derived/rqtl_gen_abh_phased.csv, rqtl_phe.csv
# Writes  data/qtl/derived/rqtl_cross_map_sweep.rds (best-supported map)
#         output/qtl/f2_map_errorprob_sweep.png
# Usage: Rscript scripts/qtl/qtl_map_errorprob_sweep.R > agent/qtl_map_errorprob_sweep.log 2>&1

suppressMessages({library(qtl); library(data.table); library(ggplot2)})
rule <- function(x) cat("\n", strrep("=", 70), "\n", x, "\n", strrep("=", 70), "\n", sep = "")

D <- "data/qtl/derived"
EPS <- c(1e-4, 1e-3, 2.5e-3, 5e-3, 1e-2, 2e-2)
NCLUST <- min(8L, parallel::detectCores())
dir.create("output/qtl", showWarnings = FALSE, recursive = TRUE)

rule("1. CROSS (phased, linkage-filtered; NO quirky)")
cr <- read.cross(format = "csvsr", dir = "",
                 genfile = file.path(D, "rqtl_gen_abh_phased.csv"),
                 phefile = file.path(D, "rqtl_phe.csv"),
                 genotypes = c("A","H","B"), na.strings = c("-","NA"),
                 crosstype = "f2", estimate.map = FALSE)
cat("markers in:", totmar(cr), " individuals:", nind(cr), "\n")

# segregation distortion, F2 expectation 0.25/0.50/0.25, per-chromosome outliers
renorm_z <- function(x) {
  z <- rep(NA_real_, length(x)); ok <- is.finite(x)
  if (sum(ok) < 5L) return(z)
  d <- ecdf(x[ok]); u <- suppressWarnings(predict(smooth.spline(x[ok], d(x[ok])), x[ok])$y)
  z[ok] <- qnorm(pmin(pmax(u, 1e-6), 1 - 1/sum(ok))); z
}
is_outlier <- function(x) { z <- renorm_z(x); !is.na(z) & z > 1.96 }
G <- pull.geno(cr); chrv <- rep(names(cr$geno), nmar(cr))
cnt <- cbind(colSums(G == 1, na.rm=TRUE), colSums(G == 2, na.rm=TRUE), colSums(G == 3, na.rm=TRUE))
tot <- rowSums(cnt); e <- tot %o% c(.25,.5,.25); chi <- rowSums((cnt - e)^2 / e)
dist_out <- logical(ncol(G))
for (ch in unique(chrv)) { i <- which(chrv == ch); dist_out[i] <- is_outlier(chi[i]) }
cat("distortion-flagged:", sum(dist_out), "\n")
cr <- pull.markers(cr, colnames(G)[!dist_out])
cat("markers used for the sweep:", totmar(cr), "\n")

rule("2. SWEEP")
res <- rbindlist(lapply(EPS, function(ep) {
  t0 <- Sys.time()
  m <- est.map(cr, error.prob = ep, map.function = "kosambi",
               maxit = 10000, tol = 1e-6, n.cluster = NCLUST)
  len <- sapply(m, function(v) max(v) - min(v))
  cat(sprintf("error.prob %-7g total %7.1f cM  (%.1f min)  per-chr: %s\n",
              ep, sum(len), as.numeric(difftime(Sys.time(), t0, units="mins")),
              paste(round(len), collapse = " ")))
  data.table(error_prob = ep, total_cM = sum(len),
             chr = names(len), cM = as.numeric(len))
}))

rule("3. SENSITIVITY -- WHERE DOES THE CURVE FLATTEN?")
tt <- unique(res[, .(error_prob, total_cM)])[order(error_prob)]
tt[, pct_change := c(NA, round(100 * diff(total_cM) / head(total_cM, -1), 1))]
print(tt)
cat("\nA flattening curve means error.prob has stopped absorbing real crossovers.\n")
cat("Expected maize F2 total: 1400-1800 cM (Ed Coe ~1781, Chen TeoNAM ~1540).\n")

p <- ggplot(tt, aes(error_prob, total_cM)) +
  geom_hline(yintercept = c(1400, 1800), linetype = 2, colour = "grey60") +
  geom_line(colour = "grey50") + geom_point(size = 2) +
  scale_x_log10() +
  labs(x = "error.prob (log scale)", y = "total map length (cM)",
       title = "Map length vs assumed genotyping-error rate",
       subtitle = "dashed band = maize F2 expectation; pick where the curve flattens, not where it hits the band") +
  theme_classic(base_size = 12)
ggsave("output/qtl/f2_map_errorprob_sweep.png", p, width = 8, height = 5, dpi = 150, bg = "white")

rule("4. PER-CHROMOSOME AT EACH error.prob")
print(dcast(res, chr ~ error_prob, value.var = "cM")[order(as.integer(chr))])

rule("5. SAVE THE MAP AT THE FLATTEST DEFENSIBLE error.prob")
# choose the smallest error.prob beyond which |%change| < 5% (curve has flattened)
flat <- tt[!is.na(pct_change) & abs(pct_change) < 5, min(error_prob)]
if (!is.finite(flat)) { flat <- min(EPS); cat("curve never flattened; using the minimum\n") }
cat("selected error.prob:", flat, "\n")
m <- est.map(cr, error.prob = flat, map.function = "kosambi",
             maxit = 10000, tol = 1e-6, n.cluster = NCLUST)
cr <- replace.map(cr, m)
saveRDS(cr, file.path(D, "rqtl_cross_map_sweep.rds"))
cat("total:", round(sum(sapply(m, function(v) max(v)-min(v))), 1), "cM over", totmar(cr), "markers\n")
cat("wrote", file.path(D, "rqtl_cross_map_sweep.rds"), "\n")
