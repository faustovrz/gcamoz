#!/usr/bin/env Rscript
# qtl_test_pileup_coords.R -- is the pileup reading the RIGHT BASE at the right place?
#
# CONTEXT: two independent routes to CML530's allele disagree on 905 of 3,115
# markers (70.95% agreement).
#   direct  (qtl_cml530_alleles_direct.R): tag alignment on CML530 + its own
#           SnpPosition -> samtools faidx one base. Internally consistent:
#           tag base == DArT ref allele in 3,730/3,730 (100%). Phase 93.3%.
#   pileup  (qtl_cml530_pileup_consensus.R): raw samtools mpileup on the
#           CML530-span-vs-B73v5 BAM at the LIFTOVER-derived v5 positions.
#
# HYPOTHESIS: the pileup is reading the wrong BASE because of a COORDINATE
# mismatch, not a base-reading error. The marker's v5 position comes from lifting
# the AGPv4 SNP coordinate; the span's v5 position comes from independently
# re-aligning CML530 sequence to v5. Those are two different derivations of
# "where is this locus in v5" and need not agree to the base, especially near
# indels. A 1 bp offset reads a neighbouring base while looking perfectly healthy.
#
# TWO TESTS
#   T1 (indirect) Does the B73 v5 base at each marker position equal the DArT REF
#      allele? DArT defined ref BY its B73 alignment, so at a correctly located
#      SNP this must hold (up to strand). Split by agree/disagree to localise.
#   T2 (direct)   Measure the offset. Take the SNP's position inside the span
#      (known from the direct route), walk the span's v5 CIGAR, and compare the
#      resulting v5 coordinate with the liftover coordinate. This is the actual
#      discrepancy in base pairs, not an inference about it.
#
# READ-ONLY.
# Usage: Rscript scripts/qtl/qtl_test_pileup_coords.R > agent/qtl_test_pileup_coords.log 2>&1

suppressMessages({library(data.table)})
rule <- function(x) cat("\n", strrep("=", 70), "\n", x, "\n", strrep("=", 70), "\n", sep = "")

D      <- "data/qtl/derived"
B73V5  <- "/Users/fvrodriguez/repos/zealtiger/data/brbseq/ref/Zm-B73-REFERENCE-NAM-5.0.fa"
BAM_V5 <- file.path(D, "CML530", "CML530_vs_B73v5.bam")
SPANS  <- file.path(D, "CML530", "CML530_spans.bed")
DIRECT <- file.path(D, "CML530_marker_alleles_direct.tsv")
PILEUP <- file.path(D, "CML530_marker_alleles_pileup.tsv")
comp   <- c(A = "T", C = "G", G = "C", T = "A")

dr <- fread(DIRECT, colClasses = list(character = "chr_v5"))
pl <- fread(PILEUP, colClasses = list(character = "chr_v5"))
cmp <- merge(dr[, .(marker, chr_v5, pos_v5, ref_allele, alt_allele,
                    cml_rname = rname, cml_refpos = ref_pos, direct = cml_class)],
             pl[, .(marker, ref_base_v5 = ref_base, pileup = cml_class)],
             by = "marker")
cmp[, agree := direct == pileup]
called <- cmp[direct %in% c("ref","alt") & pileup %in% c("ref","alt")]
cat("markers with a call from both routes:", nrow(called), "\n")
cat("agreement:", sum(called$agree), sprintf("(%.2f%%)\n", 100*mean(called$agree)))

rule("T1. DOES THE B73 v5 BASE AT THE MARKER POSITION EQUAL THE DArT REF ALLELE?")
# DArT's ref allele is B73's allele by construction, so at a correctly located SNP
# the v5 reference base must be ref_allele (or its complement if orientation flipped).
called[, v5_matches_ref := ref_base_v5 == ref_allele | ref_base_v5 == comp[ref_allele]]
called[, v5_matches_alt := ref_base_v5 == alt_allele | ref_base_v5 == comp[alt_allele]]
cat("v5 base == DArT ref allele (+/- complement):",
    sum(called$v5_matches_ref), sprintf("(%.1f%%)\n", 100*mean(called$v5_matches_ref)))
cat("v5 base == DArT alt allele (+/- complement):",
    sum(called$v5_matches_alt), sprintf("(%.1f%%)\n", 100*mean(called$v5_matches_alt)))
cat("v5 base is NEITHER allele:",
    sum(!called$v5_matches_ref & !called$v5_matches_alt),
    sprintf("(%.1f%%)\n", 100*mean(!called$v5_matches_ref & !called$v5_matches_alt)))
cat("\n(random expectation for 'matches ref' would be ~50%: 2 of 4 bases)\n")
cat("\nsplit by whether the two routes AGREE:\n")
print(called[, .(n = .N,
                 pct_v5_is_ref = round(100*mean(v5_matches_ref), 1),
                 pct_v5_is_neither = round(100*mean(!v5_matches_ref & !v5_matches_alt), 1)),
             by = agree])

rule("T2. MEASURE THE COORDINATE OFFSET DIRECTLY")
# sep="\t" is REQUIRED: the marker names are DArT AlleleIDs containing "|"
# (e.g. 7059160|F|0--28:A>G), and fread's separator auto-detection picks "|"
# over the actual tabs, producing a 6-column parse and a name-count error.
spans <- fread(SPANS, header = FALSE, sep = "\t",
               col.names = c("cml_chr","start0","end","marker"))
sam <- fread(cmd = sprintf("samtools view -F 0x904 %s | cut -f1,2,3,4,6", BAM_V5),
             header = FALSE, sep = "\t",
             col.names = c("marker","flag","v5_chr","v5_pos","cigar"))
cat("spans:", nrow(spans), " v5 alignments:", nrow(sam), "\n")
# Both files carry the FULL FASTA header as the name, not the bare AlleleID:
#   <AlleleID>|SNP=...|snppos=...|v5=chr:pos
# Strip from the unambiguous "|SNP=" delimiter. (Embedding metadata in the FASTA
# name with "|" when the AlleleID itself contains "|" was a bad design choice in
# qtl_liftover_v4_to_v5.R -- it has now caused three separate join failures.)
spans[, marker := sub("\\|SNP=.*$", "", marker)]
sam[,   marker := sub("\\|SNP=.*$", "", marker)]
x <- merge(merge(called, spans, by = "marker"), sam, by = "marker")
cat("joined for offset measurement:", nrow(x), "\n")
x[, rev_v5 := bitwAnd(flag, 16L) > 0L]
x[, span_len := end - start0]
# SNP offset inside the span, 0-based, in CML530 orientation
x[, off_in_span := cml_refpos - (start0 + 1L)]
x <- x[off_in_span >= 0 & off_in_span < span_len]
cat("with a valid in-span offset:", nrow(x), "\n")
# the span is stored in the v5 BAM in reference orientation; if it mapped reverse,
# the offset counts from the other end
x[, off_q := ifelse(rev_v5, span_len - 1L - off_in_span, off_in_span)]

q2r <- function(cigar, pos, qoff) {
  m <- regmatches(cigar, gregexpr("[0-9]+[MIDNSHP=X]", cigar))[[1]]
  len <- as.integer(sub("[MIDNSHP=X]$","",m)); op <- sub("^[0-9]+","",m)
  q <- 0L; r <- pos
  for (k in seq_along(op)) {
    o <- op[k]; L <- len[k]
    if (o %in% c("M","=","X")) { if (qoff < q + L) return(r + (qoff - q)); q <- q + L; r <- r + L }
    else if (o %in% c("I","S")) { if (qoff < q + L) return(NA_integer_); q <- q + L }
    else if (o %in% c("D","N")) { r <- r + L }
  }
  NA_integer_
}
x[, v5_from_span := mapply(q2r, cigar, v5_pos, off_q)]
x[, delta := v5_from_span - pos_v5]
cat("\nresolved:", sum(!is.na(x$delta)), "of", nrow(x), "\n")
cat("\ndelta = (v5 position implied by the span alignment) - (liftover v5 position):\n")
print(summary(x$delta))
cat("\ndelta == 0 (coordinates agree exactly):", sum(x$delta == 0, na.rm = TRUE),
    sprintf("(%.1f%%)\n", 100*mean(x$delta == 0, na.rm = TRUE)))
cat("|delta| <= 1:", sum(abs(x$delta) <= 1, na.rm = TRUE),
    sprintf("(%.1f%%)\n", 100*mean(abs(x$delta) <= 1, na.rm = TRUE)))
cat("\ndelta distribution (top 12):\n")
print(head(sort(table(x$delta), decreasing = TRUE), 12))
cat("\nsame chromosome?", sum(x$v5_chr == paste0("chr", x$chr_v5)), "of", nrow(x), "\n")

rule("T2b. OFFSET vs AGREEMENT -- the decisive cross-tab")
print(x[!is.na(delta), .(n = .N,
                         pct_delta0 = round(100*mean(delta == 0), 1),
                         median_abs_delta = median(abs(delta))),
        by = agree])

rule("VERDICT")
d0 <- mean(x$delta == 0, na.rm = TRUE)
r1 <- mean(called$v5_matches_ref)
cat(sprintf("T1: v5 base equals the DArT ref allele in %.1f%% of markers\n", 100*r1))
cat(sprintf("T2: span-implied and liftover v5 coordinates agree exactly in %.1f%%\n", 100*d0))
if (d0 < 0.9 || r1 < 0.9) {
  cat("\n=> THE PILEUP ROUTE IS MIS-LOCATED. The liftover coordinate and the span's\n")
  cat("   own v5 alignment do not point at the same base, so mpileup reads a\n")
  cat("   neighbouring position. Removing the caller fixed the ref/alt SKEW but not\n")
  cat("   this. The DIRECT route stands: it never leaves the tag's own alignment,\n")
  cat("   which is why tag base == ref allele is exactly 100%.\n")
} else {
  cat("\n=> coordinates are sound; the disagreement lies elsewhere. Investigate the\n")
  cat("   base parsing in the pileup column instead.\n")
}
