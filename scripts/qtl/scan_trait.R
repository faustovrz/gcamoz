#!/usr/bin/env Rscript
# =============================================================================
# Single-trait scanone on the current map, with permutation thresholds and
# physically-anchored support intervals.
#
#   TRAIT=TDM Rscript agent/scan_trait.R 2>&1 | tee agent/scan_TDM.log
#
# error.prob is 0.01, matching the map. Measured to be immaterial at scan time:
# 0.01 vs 2.5e-3 gave r = 0.9978 across 3,005 positions, mean dLOD +0.009
# (agent/scan_pue_errorprob.log), so no comparison arm is run here.
#
# Support intervals are reported in cM AND in CML530 Mb. The cM width is inflated with
# the map (~1,858 cM against a ~1,550 expectation); the Mb interval is not, and is what
# a candidate-gene search should use.
# =============================================================================

suppressMessages({library(qtl); library(data.table); library(ggplot2)})
rule <- function(x) cat("\n", strrep("=", 78), "\n", x, "\n", strrep("=", 78), "\n", sep = "")

D       <- "data/qtl/derived"
TRAIT   <- Sys.getenv("TRAIT", "TDM")
MODEL   <- Sys.getenv("MODEL", "normal")   # "np" for ordinal traits such as PAPP
ERROR_PROB <- 0.01
N_PERM  <- 1000
STEP    <- 1
NCLUST  <- min(8L, parallel::detectCores())
set.seed(1)

cr <- readRDS(file.path(D, "rqtl_cross_map_teonamqc.rds"))
info <- fread(file.path(D, "abh_all_marker_info.tsv"), colClasses = list(character = "chr_v5"))
stopifnot(TRAIT %in% names(cr$pheno))

rule(sprintf("1. INPUT — %s, model %s", TRAIT, MODEL))
y <- cr$pheno[[TRAIT]]
cat(sprintf("map: %d markers, %.1f cM, %d individuals\n", totmar(cr),
            sum(sapply(pull.map(cr), function(v) max(v) - min(v))), nind(cr)))
cat(sprintf("%s: n = %d non-missing of %d\n", TRAIT, sum(!is.na(y)), length(y)))
print(round(summary(y), 4))
if (MODEL == "normal") {
  sw <- shapiro.test(y[!is.na(y)])
  cat(sprintf("Shapiro-Wilk W = %.4f, p = %.3g  (normality of the raw trait)\n",
              sw$statistic, sw$p.value))
}

rule("2. SCAN + PERMUTATION")
crx <- calc.genoprob(cr, step = STEP, error.prob = ERROR_PROB, map.function = "haldane")
sc  <- scanone(crx, pheno.col = TRAIT, model = MODEL)
pm  <- scanone(crx, pheno.col = TRAIT, model = MODEL, n.perm = N_PERM,
               n.cluster = NCLUST, verbose = FALSE)
th  <- summary(pm, alpha = c(0.05, 0.10, 0.63))[, 1]
cat(sprintf("permutations: %d\n", N_PERM))
cat(sprintf("genome-wide thresholds: 5%% = %.3f   10%% = %.3f   63%% = %.3f\n",
            th[1], th[2], th[3]))

rule("3. PEAK PER CHROMOSOME")
pk <- as.data.table(sc)[, .(pos = pos[which.max(lod)], lod = max(lod)),
                        by = .(chr = as.character(chr))]
pk[, sig := fifelse(lod >= th[1], "** 0.05", fifelse(lod >= th[2], "*  0.10",
             fifelse(lod >= th[3], ".  0.63", "")))]
print(pk[order(-lod)][, .(chr, pos = round(pos, 1), lod = round(lod, 3), sig)])

rule("4. PEAKS ABOVE THE 63% THRESHOLD, WITH SUPPORT INTERVALS")
mp <- pull.map(cr)
mb_of <- function(chr, cm_lo, cm_hi) {
  v <- mp[[chr]]
  inr <- names(v)[v >= cm_lo & v <= cm_hi]
  if (!length(inr)) inr <- names(sort(abs(v - (cm_lo + cm_hi)/2)))[1]
  p <- info[marker %in% inr, ref_pos]
  c(min(p)/1e6, max(p)/1e6, length(inr))
}
hits <- pk[lod >= th[3]][order(-lod)]
if (!nrow(hits)) cat("nothing reaches even the 63% threshold.\n")
res <- list()
for (i in seq_len(nrow(hits))) {
  ch <- hits$chr[i]
  li <- lodint(sc, chr = ch, drop = 1.5, expandtomarkers = TRUE)
  mb <- mb_of(ch, min(li$pos), max(li$pos))
  # variance explained, from the single-QTL fit at the peak
  qc <- makeqtl(crx, chr = ch, pos = hits$pos[i], what = "prob")
  fq <- fitqtl(crx, pheno.col = TRAIT, qtl = qc, method = "hk", get.ests = TRUE)
  pve <- summary(fq)$result.full["Model", "%var"]
  res[[i]] <- data.table(chr = ch, peak_cM = round(hits$pos[i], 1),
                         lod = round(hits$lod[i], 3),
                         ci_lo_cM = round(min(li$pos), 1), ci_hi_cM = round(max(li$pos), 1),
                         width_cM = round(diff(range(li$pos)), 1),
                         ci_lo_Mb = round(mb[1], 2), ci_hi_Mb = round(mb[2], 2),
                         width_Mb = round(mb[2] - mb[1], 2), markers_in_ci = mb[3],
                         pct_var = round(pve, 2), sig = hits$sig[i])
}
if (length(res)) { out <- rbindlist(res); print(out) }

rule("5. TOP PEAK — nearest markers and effect")
top <- pk[which.max(lod)]
v <- mp[[top$chr]]
near <- names(sort(abs(v - top$pos)))[1:3]
print(info[marker %in% near, .(marker = substr(marker, 1, 22), cml_chr,
                               Mb = round(ref_pos/1e6, 2), v5chr = chr_v5,
                               v5_Mb = round(pos_v5/1e6, 2), maf = round(maf, 3))])
qc <- makeqtl(crx, chr = top$chr, pos = top$pos, what = "prob")
fq <- fitqtl(crx, pheno.col = TRAIT, qtl = qc, method = "hk", get.ests = TRUE)
print(summary(fq))

rule("6. OUTPUT")
dt <- as.data.table(sc); dt[, chr := as.character(chr)]
fwrite(dt, sprintf("output/qtl/scan_%s.tsv", TRAIT), sep = "\t")
dt[, chrf := factor(chr, levels = as.character(1:10))]
p <- ggplot(dt, aes(pos, lod)) +
  geom_line(linewidth = 0.4, colour = "grey20") +
  geom_hline(yintercept = th[1], linetype = "dashed", colour = "firebrick") +
  geom_hline(yintercept = th[2], linetype = "dotted", colour = "grey40") +
  facet_grid(~chrf, scales = "free_x", space = "free_x") +
  labs(x = "position (cM)", y = "LOD",
       title = sprintf("%s — scanone, error.prob %.3g, %d permutations",
                       TRAIT, ERROR_PROB, N_PERM),
       subtitle = sprintf("dashed = 5%% genome-wide (%.2f), dotted = 10%% (%.2f)",
                          th[1], th[2])) +
  theme_bw(base_size = 11) + theme(panel.spacing.x = unit(0.05, "lines"))
ggsave(sprintf("output/qtl/scan_%s.png", TRAIT), p, width = 12, height = 4, dpi = 150, bg = "white")
cat(sprintf("wrote output/qtl/scan_%s.{tsv,png}\n", TRAIT))
