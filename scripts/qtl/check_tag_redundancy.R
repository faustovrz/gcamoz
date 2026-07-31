#!/usr/bin/env Rscript
# check_tag_redundancy.R -- how many markers in the CURRENT 2,501-marker encoded set
# sit on a shared 69bp DArT tag (CloneID)?
#
# Markers on one CloneID are ONE locus: same restriction fragment, same segregation,
# same recombination information. Including several is not extra data -- it inflates
# apparent marker density and distorts the inter-marker gap distribution that
# find_quirky's data-driven threshold is computed from.
#
# This is not a QC filter. It is a statement about what a locus is.
# READ-ONLY.
# Usage: Rscript scripts/qtl/check_tag_redundancy.R > agent/check_tag_redundancy.log 2>&1

suppressMessages(library(data.table))
rule <- function(x) cat("\n", strrep("=", 70), "\n", x, "\n", strrep("=", 70), "\n", sep = "")
D <- "data/qtl/derived"

mk <- fread(file.path(D, "abh_all_marker_info.tsv"),
            colClasses = list(character = "chr_v5"))
cat("markers in the encoded set:", nrow(mk), "\n")
cat("unique CloneIDs (tags)    :", uniqueN(mk$CloneID), "\n")
cat("REDUNDANT markers         :", nrow(mk) - uniqueN(mk$CloneID), "\n\n")

tg <- mk[, .(n_snps = .N), by = CloneID]
cat("SNPs per tag:\n"); print(tg[, .N, by = n_snps][order(n_snps)])

rule("DO CO-TAG MARKERS SEGREGATE IDENTICALLY?")
# If they are truly one locus, genotypes must be near-identical. Divergence means one
# of the SNP calls on that tag is unreliable.
multi <- tg[n_snps > 1, CloneID]
cat("tags carrying >1 SNP:", length(multi), "\n")
if (length(multi)) {
  gen <- fread(file.path(D, "rqtl_gen_abh_all.csv"), sep = ",")
  setkey(gen, id)
  ids <- names(gen)[-(1:2)]
  res <- rbindlist(lapply(multi, function(cl) {
    ms <- mk[CloneID == cl, marker]
    g  <- as.matrix(gen[J(ms), ..ids])
    if (nrow(g) < 2) return(NULL)
    pr <- combn(nrow(g), 2)
    data.table(CloneID = cl, n_snps = nrow(g),
               concord = apply(pr, 2, function(p) {
                 a <- g[p[1], ]; b <- g[p[2], ]
                 ok <- !is.na(a) & !is.na(b) & a != "-" & b != "-"
                 if (sum(ok) < 30) return(NA_real_)
                 mean(a[ok] == b[ok])
               }),
               same_pos = uniqueN(mk[CloneID == cl, pos_v5]) == 1L)
  }))
  cat("\ngenotype concordance between SNPs on the same tag:\n")
  print(round(summary(res$concord), 4))
  cat("\nconcordance bands:\n")
  print(res[, .N, by = .(band = cut(concord, c(-Inf,.5,.8,.9,.95,.99,1)))][order(band)])
  cat("\npairs with concordance < 0.90:", res[concord < 0.90, .N], "of", nrow(res), "\n")
  cat("tags whose SNPs share one v5 position:", res[same_pos == TRUE, .N], "\n")
  cat("\n=> HIGH concordance confirms they are one locus and the extra copies are\n")
  cat("   redundant. LOW concordance means one SNP call on that tag is wrong, so the\n")
  cat("   tag is suspect rather than merely duplicated.\n")
}

rule("EFFECT OF COLLAPSING TO ONE SNP PER TAG")
best <- mk[order(-CallRate)][, .SD[1], by = CloneID]
cat("would retain:", nrow(best), " (drop", nrow(mk) - nrow(best), ")\n")
print(best[, .N, by = chr_v5][order(as.integer(chr_v5))])
