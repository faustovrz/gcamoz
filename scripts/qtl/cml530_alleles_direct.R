#!/usr/bin/env Rscript
# qtl_cml530_alleles_direct.R -- read CML530's allele DIRECTLY from the assembly,
# with no variant caller.
#
# WHY THIS REPLACES THE PILEUP ROUTE
# CML530 IS a parent of this F2 (per FVRZ) but was never genotyped, so its allele
# must come from the assembly. The previous pipeline
# (qtl_cml530_genotypes.sh) did: tags -> CML530 -> extract span -> map to B73 v5
# -> bcftools mpileup/call at the marker positions. Every span is ONE sequence, so
# every site has DP = 1, where `bcftools call -m` has no power and defaults toward
# the reference. Result: 2,086 "ref" vs 1,033 "alt" -- a 2:1 skew that cannot be
# real (B73's allele is arbitrary with respect to which parent carries it), and
# polarity that was only 53.4% phase-consistent, i.e. random
# (agent/qtl_diagnose_phase.log).
#
# An assembly has ONE deterministic base at one position. There is nothing to
# infer. So: take the tag alignment to CML530, walk the CIGAR to the SNP offset,
# and read the CML530 reference base at that position with samtools faidx.
#
# Everything stays in REFERENCE orientation, so no complementing is needed for the
# ref/alt decision: the aligned tag base and the CML530 base are directly
# comparable. TrimmedSequenceRef carries the DArT REFERENCE allele, so
#   CML530 base == aligned tag base  ->  CML530 carries the ref allele
#   CML530 base != aligned tag base  ->  CML530 carries a non-ref allele
# and we additionally check that the non-ref base equals the DArT alt allele
# (strand-corrected) rather than a third base.
#
# Reads   data/qtl/derived/CML530/tags_vs_CML530.bam
#         data/qtl/derived/markers_v5.tsv
#         zealhmm .../Zm-CML530-REFERENCE-HiLo-1.0.fa.gz  (BGZF, faidx in place)
# Writes  data/qtl/derived/CML530_marker_alleles_direct.tsv
# Usage: Rscript scripts/qtl/qtl_cml530_alleles_direct.R > agent/qtl_cml530_alleles_direct.log 2>&1

suppressMessages({library(data.table)})
rule <- function(x) cat("\n", strrep("=", 70), "\n", x, "\n", strrep("=", 70), "\n", sep = "")

D      <- "data/qtl/derived"
BAM    <- file.path(D, "CML530", "tags_vs_CML530.bam")
ROSTER <- file.path(D, "markers_v5.tsv")
VCF    <- file.path(D, "CML530", "CML530_at_markers.vcf")
CML_GZ <- "/Users/fvrodriguez/repos/zealhmm/data/ref/genomes/Zm-CML530-REFERENCE-HiLo-1.0/Zm-CML530-REFERENCE-HiLo-1.0.fa.gz"
OUT    <- file.path(D, "CML530_marker_alleles_direct.tsv")

rule("1. CONFIRM THE DEPTH-1 CALLING BIAS")
dp <- system2("bcftools", c("query", "-f", "'[%DP]\\n'", VCF), stdout = TRUE, stderr = FALSE)
dp <- suppressWarnings(as.integer(gsub("'", "", dp)))
cat("pileup depth at called sites:\n"); print(table(dp, useNA = "ifany"))
cat(sprintf("\n=> %.1f%% of sites are DEPTH 1 (a few are 2-4 where tags overlap).\n",
            100 * mean(dp == 1, na.rm = TRUE)))
cat("   bcftools call -m has essentially no discriminating power at DP=1 and\n")
cat("   defaults toward REF, producing the 2:1 ref/alt skew.\n")

rule("2. PARSE THE TAG ALIGNMENTS")
roster <- fread(ROSTER, colClasses = list(character = "chr_v5"))
roster[, ref_allele := sub(".*:([ACGT])>([ACGT]).*", "\\1", SNP)]
roster[, alt_allele := sub(".*:([ACGT])>([ACGT]).*", "\\2", SNP)]
sam <- fread(cmd = sprintf("samtools view -F 0x904 %s | cut -f1,2,3,4,6,10", BAM),
             header = FALSE, sep = "\t",
             col.names = c("marker", "flag", "rname", "pos", "cigar", "seq"))
cat("primary aligned tags:", nrow(sam), "\n")
# The tag FASTA headers written by qtl_liftover_v4_to_v5.R are
#   <AlleleID>|SNP=...|snppos=...|v5=chr:pos
# with NO whitespace, so samtools takes the whole string as the read name. The
# AlleleID itself contains "|" (e.g. "7059160|F|0--28:A>G"), so splitting on "|"
# is ambiguous -- strip from the unambiguous "|SNP=" delimiter instead.
cat("raw read name example:", sam$marker[1], "\n")
sam[, marker := sub("\\|SNP=.*$", "", marker)]
cat("parsed AlleleID     :", sam$marker[1], "\n")
sam <- merge(sam, roster[, .(marker, SnpPosition, ref_allele, alt_allele,
                             chr_v5, pos_v5, framework, qc_pass)],
             by = "marker")
cat("joined to roster    :", nrow(sam), "\n")
sam[, rev := bitwAnd(flag, 16L) > 0L]
cat("reverse-strand alignments:", sum(sam$rev), "\n")

rule("3. WALK THE CIGAR: SNP OFFSET -> REFERENCE POSITION")
# BAM SEQ is always stored in reference orientation. SnpPosition is a 0-based
# offset in the ORIGINAL tag, so on a minus-strand alignment the SNP sits at
# (tag_len - 1 - SnpPosition) within SEQ.
sam[, taglen := nchar(seq)]
sam[, snp_in_seq := ifelse(rev, taglen - 1L - as.integer(SnpPosition),
                                as.integer(SnpPosition))]

# For each read, map a query offset to a reference position using the CIGAR.
q2r <- function(cigar, pos, qoff) {
  m <- regmatches(cigar, gregexpr("[0-9]+[MIDNSHP=X]", cigar))[[1]]
  len <- as.integer(sub("[MIDNSHP=X]$", "", m)); op <- sub("^[0-9]+", "", m)
  q <- 0L; r <- pos                      # q: 0-based in SEQ, r: 1-based ref
  for (k in seq_along(op)) {
    o <- op[k]; L <- len[k]
    if (o %in% c("M", "=", "X")) {
      if (qoff < q + L) return(r + (qoff - q))
      q <- q + L; r <- r + L
    } else if (o %in% c("I", "S")) {
      if (qoff < q + L) return(NA_integer_)   # SNP inside an insertion/clip
      q <- q + L
    } else if (o %in% c("D", "N")) {
      r <- r + L
    }
  }
  NA_integer_
}
sam[, ref_pos := mapply(q2r, cigar, pos, snp_in_seq)]
cat("SNP offset resolved to a reference base:", sum(!is.na(sam$ref_pos)),
    "of", nrow(sam), "\n")
cat("unresolved (SNP fell in an insertion or soft clip):", sum(is.na(sam$ref_pos)), "\n")
# the aligned tag's own base at that offset, in reference orientation
sam[, tag_base := substr(seq, snp_in_seq + 1L, snp_in_seq + 1L)]
ok <- sam[!is.na(ref_pos) & tag_base %in% c("A", "C", "G", "T")]
cat("usable alignments:", nrow(ok), "\n")

rule("4. READ THE CML530 BASE (samtools faidx, one base each)")
reg <- sprintf("%s:%d-%d", ok$rname, ok$ref_pos, ok$ref_pos)
rf <- tempfile(); writeLines(reg, rf)
fa <- system2("samtools", c("faidx", "-r", rf, CML_GZ), stdout = TRUE)
# faidx emits ">region" then the single base, in the order requested
hdr <- grepl("^>", fa)
bases <- toupper(fa[!hdr])
stopifnot(sum(hdr) == nrow(ok), length(bases) == nrow(ok))
ok[, cml_base := bases]
unlink(rf)
cat("bases read:", nrow(ok), " composition:\n"); print(table(ok$cml_base))

rule("5. CLASSIFY")
comp <- c(A = "T", C = "G", G = "C", T = "A")
# alt allele expressed in reference orientation
ok[, alt_refor := ifelse(rev, comp[alt_allele], alt_allele)]
ok[, ref_refor := ifelse(rev, comp[ref_allele], ref_allele)]
# sanity: the aligned tag base SHOULD be the ref allele in reference orientation
cat("tag base == ref allele (strand-corrected):",
    sum(ok$tag_base == ok$ref_refor), "of", nrow(ok),
    sprintf("(%.1f%%)\n", 100 * mean(ok$tag_base == ok$ref_refor)))
ok[, cml_class := fifelse(cml_base == ref_refor, "ref",
                   fifelse(cml_base == alt_refor, "alt", "third-allele"))]
print(ok[, .N, by = cml_class][order(-N)])
cat("\nref:alt ratio =", round(ok[cml_class == "ref", .N] / ok[cml_class == "alt", .N], 3),
    " (a true parent should give ~1.0; the pileup route gave 2.02)\n")
cat("\nframework subset:\n")
print(ok[framework == TRUE, .N, by = cml_class][order(-N)])

rule("6. WRITE")
res <- merge(roster[, .(marker, CloneID, chr_v4, pos_v4, chr_v5, pos_v5, SNP,
                        SnpPosition, ref_allele, alt_allele, qc_pass, framework)],
             ok[, .(marker, rname, ref_pos, rev, tag_base, cml_base,
                    ref_refor, alt_refor, cml_class)],
             by = "marker", all.x = TRUE)
res[is.na(cml_class), cml_class := "no-hit"]
setorder(res, chr_v5, pos_v5)
fwrite(res, OUT, sep = "\t")
cat("wrote", OUT, " rows:", nrow(res), "\n")
print(res[, .N, by = cml_class][order(-N)])
cat("\nNEXT: re-run qtl_make_abh.R against this file and re-test phase\n")
cat("consistency. If CML530 is the parent and this reading is correct, markers\n")
cat("polarized from it should now be ~100% in phase, not 53%.\n")
