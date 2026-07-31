#!/usr/bin/env Rscript
# qtl_marker_accounting.R -- exactly what marker set is the map built on, and what
# removed the rest. Every count read from the artefacts, not recited from memory.
# READ-ONLY.
# Usage: Rscript scripts/qtl/qtl_marker_accounting.R > agent/qtl_marker_accounting.log 2>&1

suppressMessages({library(qtl); library(data.table)})
rule <- function(x) cat("\n", strrep("=", 70), "\n", x, "\n", strrep("=", 70), "\n", sep = "")
D <- "data/qtl/derived"
n <- function(f, ...) if (file.exists(f)) nrow(fread(f, ...)) else NA_integer_

rule("1. THE FUNNEL")
qcall <- fread(file.path(D, "markers_v5_qc.csv"))
ros   <- fread(file.path(D, "markers_v5.tsv"), colClasses = list(character="chr_v5"))
abh   <- fread(file.path(D, "abh_marker_info.tsv"), colClasses = list(character="chr_v5"))
gen   <- fread(file.path(D, "rqtl_gen_abh.csv"), sep = ",", select = 1:2)
phz   <- if (file.exists(file.path(D, "rqtl_gen_abh_phased.csv")))
           fread(file.path(D, "rqtl_gen_abh_phased.csv"), sep = ",", select = 1:2) else NULL
crmap <- if (file.exists(file.path(D, "rqtl_cross_map.rds")))
           readRDS(file.path(D, "rqtl_cross_map.rds")) else NULL

fun <- data.table(
  step = c("SNPs in the DArT report",
           "v4-anchored to chr 1-10",
           "lifted 1:1 to v5, same chromosome",
           "qc_pass (CallRate>=.90, RepAvg>=.95, MAF>=.15, 1:2:1 p>.01)",
           "framework (one SNP per 69bp tag)",
           "ABH-encoded genotype file",
           "phased file  [MY LINK_MIN filter applied here]",
           "in the current fitted map"),
  markers = c(nrow(qcall),
              qcall[, sum(chr_v4 %in% as.character(1:10))],
              nrow(ros),
              ros[, sum(qc_pass)],
              ros[, sum(framework)],
              nrow(gen),
              if (!is.null(phz)) nrow(phz) else NA_integer_,
              if (!is.null(crmap)) totmar(crmap) else NA_integer_))
fun[, removed := c(NA, -diff(markers))]
print(fun)

rule("2. *** WHAT THE 'AIRMINE FILTERS ONLY' RUN WOULD ACTUALLY USE ***")
cat("The phased genotype file already has MY invented LINK_MIN filter baked in:\n")
cat(sprintf("  ABH-encoded : %d markers\n  phased file : %d markers  (%d removed by LINK_MIN |cor| < 0.80)\n",
            nrow(gen), if (!is.null(phz)) nrow(phz) else NA,
            nrow(gen) - (if (!is.null(phz)) nrow(phz) else NA)))
cat("\nSo running 'their filters only' against the phased file is NOT their filters\n")
cat("only -- it inherits mine. To test their pipeline honestly the phasing step must\n")
cat("re-run WITHOUT the marker drop, keeping all 1,384, and let distortion +\n")
cat("find_quirky be the only filters.\n")

rule("3. WHAT REMOVED MARKERS AFTER THE FRAMEWORK STAGE")
if (file.exists(file.path(D, "map_qc.tsv"))) {
  mq <- fread(file.path(D, "map_qc.tsv"), colClasses = list(character="chr_v5"))
  print(mq[, .(n = .N), by = .(dropped_distortion, dropped_quirky, in_map)][order(-n)])
  cat("\nper chromosome, markers surviving into the map:\n")
  print(mq[, .(framework = .N, in_map = sum(in_map),
               dist = sum(dropped_distortion), quirky = sum(dropped_quirky)),
           by = chr_v5][order(as.integer(chr_v5))])
}

rule("4. COMPOSITION OF THE CURRENT MAP")
if (!is.null(crmap)) {
  inmap <- unlist(lapply(pull.map(crmap), names))
  a <- abh[marker %in% inmap]
  cat("markers in map:", length(inmap), "\n\npolarity source:\n")
  print(a[, .N, by = polarity_source])
  cat("\nCML530 allele class:\n"); print(a[, .N, by = cml_class][order(-N)])
  cat("\nmarker quality of those in the map:\n")
  print(a[, .(median_CallRate = round(median(CallRate),3),
              median_RepAvg   = round(median(RepAvg),3),
              median_maf      = round(median(maf),3),
              median_p121     = signif(median(p121),3))])
  cat("\nsame, for framework markers NOT in the map:\n")
  b <- abh[!marker %in% inmap]
  print(b[, .(n = .N,
              median_CallRate = round(median(CallRate),3),
              median_RepAvg   = round(median(RepAvg),3),
              median_maf      = round(median(maf),3),
              median_p121     = signif(median(p121),3))])
  cat("\n=> if the excluded markers look as good as the included ones on CallRate,\n")
  cat("   RepAvg, MAF and 1:2:1 fit, they were dropped for MAP-GEOMETRY reasons\n")
  cat("   (find_quirky), not for data quality.\n")
}

rule("5. FILTER PROVENANCE")
cat("THEIRS (airmine / TeoNAM):\n")
cat("  * segregation distortion, per-chromosome renorm_z outliers (z>1.96)\n")
cat("  * find_quirky: data-driven fine gap threshold + island rule (2 cM, max 20)\n")
cat("MINE (invented this session):\n")
cat("  * qc_pass marker filter (CallRate/RepAvg/MAF/1:2:1) -- upstream, at ABH build\n")
cat("  * framework reduction to one SNP per tag\n")
cat("  * LINK_MIN local-linkage drop, |cor| < 0.80 with +/-5 neighbours\n")
cat("  * phase correction (windowed sign flip) -- NOT a filter, but it changes coding\n")
cat("\nNOTE the qc_pass + framework steps are also mine and they are the LARGEST\n")
cat("reduction of all: 3,828 -> 1,384. Any comparison of filter regimes has to say\n")
cat("whether it starts from 3,828 or from 1,384.\n")
