#!/usr/bin/env Rscript
# =============================================================================
# scanone for every F2 phenotype, with per-trait permutation thresholds.
#
#   Rscript agent/scan_all_traits.R 2>&1 | tee agent/scan_all_traits.log
#
# error.prob = 0.01 throughout, matching the map (measured immaterial at scan time:
# r = 0.9978 vs 2.5e-3, agent/scan_pue_errorprob.log).
#
# PAPP is an ordinal 1-4 score, so it is scanned with model = "np" (Kruskal-Wallis);
# everything else uses the normal model. Normality is REPORTED per trait because
# several of these are strongly skewed and that inflates permutation thresholds.
#
# Multiple traits are scanned, so a per-trait 5% threshold is NOT a 5% experiment-wide
# rate. The trait-corrected threshold is reported alongside.
# =============================================================================

suppressMessages({library(qtl); library(data.table)})
rule <- function(x) cat("\n", strrep("=", 78), "\n", x, "\n", strrep("=", 78), "\n", sep = "")

D <- "data/qtl/derived"
ERROR_PROB <- 0.01; N_PERM <- 1000; STEP <- 1
NCLUST <- min(8L, parallel::detectCores())
set.seed(1)

cr <- readRDS(file.path(D, "rqtl_cross_map_teonamqc.rds"))
info <- fread(file.path(D, "abh_all_marker_info.tsv"), colClasses = list(character = "chr_v5"))
crx <- calc.genoprob(cr, step = STEP, error.prob = ERROR_PROB, map.function = "haldane")
mp <- pull.map(cr)

TRAITS <- setdiff(names(cr$pheno), "id")
NP <- "PAPP"                      # ordinal
cat(sprintf("map: %d markers, %.1f cM | traits: %s\n", totmar(cr),
            sum(sapply(mp, function(v) max(v) - min(v))), paste(TRAITS, collapse = " ")))

mb_at <- function(chr, cm) {
  v <- mp[[chr]]
  as.numeric(approx(v, info$ref_pos[match(names(v), info$marker)], xout = cm, rule = 2)$y) / 1e6
}

rows <- list(); prof <- list()
for (tr in TRAITS) {
  model <- if (tr %in% NP) "np" else "normal"
  y <- cr$pheno[[tr]]
  n_ok <- sum(!is.na(y))
  if (n_ok < 50 || length(unique(y[!is.na(y)])) < 3) {
    cat(sprintf("%-6s SKIPPED (n = %d, %d distinct values)\n", tr, n_ok,
                length(unique(y[!is.na(y)]))))
    next
  }
  sw <- if (model == "normal") shapiro.test(y[!is.na(y)])$p.value else NA_real_
  sc <- scanone(crx, pheno.col = tr, model = model)
  pm <- scanone(crx, pheno.col = tr, model = model, n.perm = N_PERM,
                n.cluster = NCLUST, verbose = FALSE)
  th <- summary(pm, alpha = c(0.05, 0.10))[, 1]
  i  <- which.max(sc$lod)
  pk_chr <- as.character(sc$chr[i]); pk_cm <- sc$pos[i]; pk_lod <- sc$lod[i]

  # contiguous 1.5-LOD interval around the peak
  scc <- as.data.table(sc)[as.character(chr) == pk_chr][order(pos)]
  j <- which.max(scc$lod); cut <- pk_lod - 1.5
  lo <- j; while (lo > 1 && scc$lod[lo - 1] >= cut) lo <- lo - 1
  hi <- j; while (hi < nrow(scc) && scc$lod[hi + 1] >= cut) hi <- hi + 1

  rows[[length(rows) + 1L]] <- data.table(
    trait = tr, model = model, n = n_ok, shapiro_p = signif(sw, 3),
    peak_chr = pk_chr, peak_cM = round(pk_cm, 1),
    peak_Mb = round(mb_at(pk_chr, pk_cm), 2), lod = round(pk_lod, 3),
    thr05 = round(th[1], 3), thr10 = round(th[2], 3),
    sig05 = pk_lod >= th[1], sig10 = pk_lod >= th[2],
    ci_lo_cM = round(scc$pos[lo], 1), ci_hi_cM = round(scc$pos[hi], 1),
    ci_lo_Mb = round(mb_at(pk_chr, scc$pos[lo]), 2),
    ci_hi_Mb = round(mb_at(pk_chr, scc$pos[hi]), 2))
  d <- as.data.table(sc); d[, `:=`(trait = tr, chr = as.character(chr))]
  prof[[length(prof) + 1L]] <- d
  r <- rows[[length(rows)]]
  cat(sprintf("%-6s %-6s n=%3d | peak chr%-2s %6.1f cM (%7.2f Mb) LOD %5.3f | thr05 %5.3f | %s\n",
              tr, model, n_ok, r$peak_chr, r$peak_cM, r$peak_Mb, r$lod, r$thr05,
              ifelse(r$sig05, "SIGNIFICANT", ifelse(r$sig10, "suggestive", "-"))))
}
res <- rbindlist(rows); pr <- rbindlist(prof)

rule("SUMMARY — all traits, ranked by LOD")
print(res[order(-lod), .(trait, model, n, peak_chr, peak_cM, peak_Mb, lod,
                         thr05, sig05, sig10)])

rule("TRAIT-CORRECTED SIGNIFICANCE")
k <- nrow(res)
cat(sprintf("%d traits scanned. A per-trait 5%% threshold gives an experiment-wide\n", k))
cat(sprintf("false-positive rate of 1-0.95^%d = %.1f%% if the traits were independent.\n",
            k, 100 * (1 - 0.95^k)))
cat("They are NOT independent -- several are arithmetic functions of others\n")
cat("(GY, HI, PUE, RSR, RYE, TDM), so the true rate is between 5% and that figure.\n")
cat(sprintf("\nBonferroni-corrected per-trait alpha: %.4f\n", 0.05 / k))
res[, sig05_bonf := lod >= thr05 + log10(k)]   # crude: shift on the LOD scale
print(res[sig05 == TRUE | sig10 == TRUE,
          .(trait, peak_chr, peak_Mb, lod, thr05, sig05, sig10)])

rule("SUPPORT INTERVALS FOR ANYTHING SUGGESTIVE OR BETTER")
hits <- res[sig10 == TRUE]
if (nrow(hits)) print(hits[, .(trait, peak_chr, peak_cM, lod,
                               ci_lo_cM, ci_hi_cM, ci_lo_Mb, ci_hi_Mb,
                               width_Mb = round(ci_hi_Mb - ci_lo_Mb, 1))])
if (!nrow(hits)) cat("nothing reaches even the 10% genome-wide threshold.\n")

rule("OUTPUT")
fwrite(res, "output/qtl/scan_all_summary.tsv", sep = "\t")
fwrite(pr,  "output/qtl/scan_all_profiles.tsv", sep = "\t")
cat("wrote output/qtl/scan_all_{summary,profiles}.tsv\n")

# PUBLIC. data/qtl/ is gitignored (raw genotypes), so the scan RESULTS -- no individual
# data, just LOD by position and the per-trait permutation thresholds -- go to a tracked
# path alongside data/f2_genetic_map.tsv.
pub <- copy(res)
setnames(pub, "thr05", "perm_thr_0.05"); setnames(pub, "thr10", "perm_thr_0.10")
fwrite(pub, "data/f2_qtl_scan_summary.tsv", sep = "\t")
prp <- copy(pr)[, .(trait, chr, pos_cM = round(pos, 2), lod = round(lod, 4))]
fwrite(prp, "data/f2_qtl_scan_profiles.tsv", sep = "\t")
cat(sprintf("wrote data/f2_qtl_scan_summary.tsv (%d traits) and\n", nrow(pub)))
cat(sprintf("      data/f2_qtl_scan_profiles.tsv (%d positions x trait)\n", nrow(prp)))
