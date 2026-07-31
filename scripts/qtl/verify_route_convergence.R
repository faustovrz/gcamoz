#!/usr/bin/env Rscript
# verify_route_convergence.R -- the two routes to CML530's allele MUST agree.
# If they do not, a coordinate mapping is wrong, marker positions are wrong, and the
# map is wrong. This is a MAP-VALIDITY test, not a bookkeeping curiosity.
#
# WHY THEY MUST CONVERGE (my earlier framing was wrong):
#   direct route  -- read CML530's base at the tag's SnpPosition via the tag's own
#                    alignment to CML530, compare to DArT's ref/alt
#   pileup route  -- map the CML530 span to B73 v5, pileup at the marker's LIFTED v5
#                    position, read the observed base
# Both read THE SAME BASE AT THE SAME LOCUS. I previously classified the pileup result
# as "does CML530 match B73?", which is a different question and produced only 71%
# agreement. That was my comparison being wrong, not the routes measuring different
# things. Compare the pileup's OBSERVED BASE to DArT's ref/alt -- as the direct route
# does -- and they must converge.
#
# Residual disagreement localises a real fault:
#   * liftover put the marker at the wrong v5 base, OR
#   * the span's independent v5 alignment is offset, OR
#   * the tag's CML530 alignment is paralogous
# Any of those means the physical order feeding est.map is wrong somewhere.
#
# READ-ONLY.
# Usage: Rscript scripts/qtl/verify_route_convergence.R > agent/verify_route_convergence.log 2>&1

suppressMessages(library(data.table))
rule <- function(x) cat("\n", strrep("=", 70), "\n", x, "\n", strrep("=", 70), "\n", sep = "")

D      <- "data/qtl/derived"
B73V5  <- "/Users/fvrodriguez/repos/zealtiger/data/brbseq/ref/Zm-B73-REFERENCE-NAM-5.0.fa"
BAM_V5 <- file.path(D, "CML530", "CML530_vs_B73v5.bam")
BED    <- file.path(D, "markers_v5.bed")
DIRECT <- file.path(D, "CML530_marker_alleles_direct.tsv")
comp   <- c(A = "T", C = "G", G = "C", T = "A")

rule("1. PILEUP: OBSERVED CML530 BASE AT THE LIFTED v5 POSITION")
stopifnot(file.exists(BAM_V5), file.exists(BED))
pu <- fread(cmd = sprintf("samtools mpileup -B -Q0 -q0 -f %s -l %s %s 2>/dev/null",
                          B73V5, BED, BAM_V5),
            header = FALSE, sep = "\t", fill = TRUE,
            col.names = c("chrom","pos","ref_base","depth","bases","quals"))
cat("pileup sites:", nrow(pu), "\n")
# '.' / ',' mean the span matches the B73 reference base; anything else is the
# observed base. Both are in B73 v5 orientation.
pu[, b := toupper(gsub("[^ACGTacGT.,]", "", bases))]
pu[, cml_base_v5 := fifelse(grepl("^[.,]+$", b) | b == "",
                            toupper(ref_base), toupper(gsub("[.,]", "", b)))]
pu[, cml_base_v5 := substr(cml_base_v5, 1, 1)]
pu <- pu[cml_base_v5 %in% c("A","C","G","T")]
cat("usable:", nrow(pu), "\n")

rule("2. CLASSIFY THE PILEUP BASE AGAINST DArT's ref/alt  (both strands threaded)")
dr <- fread(DIRECT, colClasses = list(character = "chr_v5"))
dr[, chrom := paste0("chr", chr_v5)]
# TWO strand flips compose between the DArT alleles and the v5 pileup base:
#   DArT alleles are in TAG orientation
#     -> tag aligned to CML530        strand A = `rev`     (tags_vs_CML530.bam)
#     -> span extracted in CML530 orientation
#     -> span aligned to B73 v5       strand B = `rev_v5`  (CML530_vs_B73v5.bam)
#     -> mpileup reports bases in v5 REFERENCE orientation
# Net orientation is therefore rev XOR rev_v5. Complement the DArT alleles when TRUE.
#
# The previous version accepted EITHER orientation, which (a) made A/T and C/G SNPs
# unresolvable -- for those the complement of one allele IS the other -- and (b) let a
# WRONG base pass whenever it equalled the complement of the right one, so that 98.32%
# was an optimistic bound rather than a tight one.
sam <- fread(cmd = sprintf("samtools view -F 0x904 %s | cut -f1,2", BAM_V5),
             header = FALSE, sep = "\t", col.names = c("marker","flag"))
sam[, marker := sub("\\|SNP=.*$", "", marker)]     # headers carry appended metadata
sam[, rev_v5 := bitwAnd(flag, 16L) > 0L]
cat("v5 alignments read:", nrow(sam), " reverse-strand:", sum(sam$rev_v5), "\n")

x <- merge(dr[, .(marker, chrom, pos_v5, ref_allele, alt_allele,
                  rev_cml = rev, direct = cml_class)],
           pu[, .(chrom, pos_v5 = pos, cml_base_v5, depth)],
           by = c("chrom","pos_v5"))
x <- merge(x, sam[, .(marker, rev_v5)], by = "marker")
cat("markers with both routes and a v5 strand:", nrow(x), "\n")
x[, net_rev := xor(rev_cml, rev_v5)]
cat("net orientation reversed:", sum(x$net_rev), "of", nrow(x), "\n")
x[, ref_or := fifelse(net_rev, comp[ref_allele], ref_allele)]
x[, alt_or := fifelse(net_rev, comp[alt_allele], alt_allele)]
x[, pileup := fifelse(cml_base_v5 == ref_or, "ref",
               fifelse(cml_base_v5 == alt_or, "alt", "neither"))]
print(x[, .N, by = pileup][order(-N)])
cat("\nNo 'ambiguous' class now: with the net strand known, each DArT allele has ONE\n")
cat("orientation, so A/T and C/G SNPs resolve like any other.\n")

rule("3. CONVERGENCE")
cmp <- x[pileup %in% c("ref","alt") & direct %in% c("ref","alt")]
cat("comparable markers:", nrow(cmp), "\n")
agree <- cmp[pileup == direct, .N]
cat(sprintf("AGREE: %d of %d (%.2f%%)\n", agree, nrow(cmp), 100*agree/nrow(cmp)))
cat("\nconfusion matrix (rows pileup, cols direct):\n")
print(table(pileup = cmp$pileup, direct = cmp$direct))
cat("\nreference points: comparing to B73 gave 70.95%; accepting either orientation gave 98.32%\n")

rule("4. WHERE THEY STILL DISAGREE")
bad <- cmp[pileup != direct]
cat("disagreements:", nrow(bad), "\n")
if (nrow(bad)) {
  cat("\nper chromosome:\n"); print(bad[, .N, by = chrom][order(-N)])
  cat("\nby pileup depth:\n"); print(bad[, .N, by = depth][order(depth)])
  cat("\nby CML530-alignment strand (direct route):\n"); print(bad[, .N, by = rev_cml])
  cat("\nfirst 15:\n")
  print(head(bad[, .(marker, chrom, pos_v5, ref_allele, alt_allele,
                     cml_base_v5, pileup, direct)], 15))
  # PERSIST. These markers have unreliable physical positions, so they are placed at
  # the wrong point in the map order and will show poor linkage to whatever they sit
  # next to. Previously they were printed and discarded, which made that untestable.
  fwrite(bad[, .(marker, chrom, pos_v5, ref_allele, alt_allele, cml_base_v5,
                 depth, rev_cml, rev_v5, pileup, direct)],
         file.path(D, "route_disagreements.tsv"), sep = "\t")
  cat(sprintf("\nwrote %s (%d markers)\n",
              file.path(D, "route_disagreements.tsv"), nrow(bad)))
}

rule("VERDICT")
pct <- 100*agree/nrow(cmp)
if (pct >= 99) {
  cat(sprintf("CONVERGED (%.2f%%). Two independent coordinate paths -- the AGPv4->v5\n", pct))
  cat("liftover, and the span's own re-alignment to v5 -- put the marker on the same\n")
  cat("base and read the same allele. The physical positions feeding est.map are sound.\n")
} else {
  cat(sprintf("NOT CONVERGED (%.2f%%). A coordinate mapping is wrong for the %d\n",
              pct, nrow(bad)))
  cat("disagreeing markers, so their physical positions are unreliable and any map\n")
  cat("built on them is suspect. Section 4 localises them; drop them or resolve the\n")
  cat("coordinate discrepancy before mapping.\n")
}
