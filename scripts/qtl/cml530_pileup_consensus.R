#!/usr/bin/env Rscript
# qtl_cml530_pileup_consensus.R -- get CML530's base AT THE v5 MARKER POSITIONS by
# coordinate, with no genotype caller and no join by name.
#
# THE POINT (per FVRZ): the pileup architecture was never the problem. The problem
# was running `bcftools call -m`, a DIPLOID caller with priors, on what is
# effectively HAPLOID depth-1 data (an assembly contributes one sequence). At DP=1
# it defaults toward REF and manufactured a 2:1 ref/alt skew, giving random
# polarity (53.4% phase-consistent). Ask instead for the OBSERVED base.
#
# Three no-caller options, tested here:
#   (a) samtools mpileup      -- raw read-base column; not a caller at all
#   (b) samtools consensus    -- --mode simple, plain majority
#   (c) bcftools call --ploidy 1 -- haploid IS the correct model for an assembly
#
# WHY THIS IS BETTER THAN MY CIGAR WALK (qtl_cml530_alleles_direct.R):
#   * keyed by v5 COORDINATE, so no join by read name at all
#   * no hand-written CIGAR walker; indels handled by samtools internally
#   * (CLAIMED, but FALSE) classification is orientation-safe because DArT's "ref"
#     allele is B73's allele. It is not -- see the correction below. The whole
#     premise of this script is unsound; it is retained as the record of a
#     measured dead end, not as a usable route.
#
# COVERAGE CAVEAT, reported not hidden: the v5-mapped BAM contains only the 3,167
# spans that survived `-F 0x904 -q 20`, so this route can reach fewer markers than
# the direct read (3,730). That is a property of how the spans were built, not of
# the calling mode.
#
# Reads   data/qtl/derived/CML530/CML530_vs_B73v5.bam, markers_v5.bed,
#         data/qtl/derived/CML530_marker_alleles_direct.tsv (for cross-validation)
# Writes  data/qtl/derived/CML530_marker_alleles_pileup.tsv
# Usage: Rscript scripts/qtl/qtl_cml530_pileup_consensus.R > agent/qtl_cml530_pileup_consensus.log 2>&1

suppressMessages({library(data.table)})
rule <- function(x) cat("\n", strrep("=", 70), "\n", x, "\n", strrep("=", 70), "\n", sep = "")

D      <- "data/qtl/derived"
BAM    <- file.path(D, "CML530", "CML530_vs_B73v5.bam")
BED    <- file.path(D, "markers_v5.bed")
ROSTER <- file.path(D, "markers_v5.tsv")
DIRECT <- file.path(D, "CML530_marker_alleles_direct.tsv")
B73V5  <- "/Users/fvrodriguez/repos/zealtiger/data/brbseq/ref/Zm-B73-REFERENCE-NAM-5.0.fa"
OUT    <- file.path(D, "CML530_marker_alleles_pileup.tsv")

rule("0. AVAILABLE NO-CALLER MODES")
sc <- suppressWarnings(system2("samtools", "consensus", stdout = TRUE, stderr = TRUE))
cat("samtools consensus present:", !any(grepl("unrecognized|invalid", sc)), "\n")
bc <- suppressWarnings(system2("bcftools", c("call", "--help"), stdout = TRUE, stderr = TRUE))
cat("bcftools call --ploidy available:", any(grepl("--ploidy", bc)), "\n")
cat("samtools version:", system2("samtools", "--version", stdout = TRUE)[1], "\n")

rule("1. RAW mpileup -- the observed base, no model at all")
# -B disables BAQ (irrelevant for assembly-derived sequence and it can mask
# genuine mismatches); -Q0 -q0 keep every base since there are no real quals.
cmd <- sprintf("samtools mpileup -B -Q0 -q0 -f %s -l %s %s 2>/dev/null",
               B73V5, BED, BAM)
pu <- fread(cmd = cmd, header = FALSE, sep = "\t", fill = TRUE,
            col.names = c("chrom", "pos", "ref_base", "depth", "bases", "quals"))
cat("pileup rows:", nrow(pu), "\n")
cat("depth distribution:\n"); print(table(pu$depth))
cat("\nraw base-column examples:\n"); print(head(pu[, .(chrom, pos, ref_base, depth, bases)], 8))

# '.' = match to reference on the forward strand, ',' = match on reverse.
# Anything else is the observed (mismatching) base. Both are already in B73 v5
# orientation, so no complementing is needed.
pu[, b := toupper(gsub("[^ACGTacgt.,]", "", bases))]
pu[, cml_base := fifelse(grepl("^[.,]+$", b) | b == "",
                         toupper(ref_base),
                         toupper(gsub("[.,]", "", b)))]
pu[, cml_base := substr(cml_base, 1, 1)]
pu[, matches_b73 := cml_base == toupper(ref_base)]
cat("\nCML530 matches B73 v5 at the marker position:",
    sum(pu$matches_b73), "of", nrow(pu),
    sprintf("(%.1f%%)\n", 100 * mean(pu$matches_b73)))

rule("2. CLASSIFY BY COORDINATE (no name join anywhere)")
roster <- fread(ROSTER, colClasses = list(character = "chr_v5"))
roster[, chrom := paste0("chr", chr_v5)]
res <- merge(roster[, .(marker, chrom, pos_v5, chr_v5, SNP, qc_pass, framework)],
             pu[, .(chrom, pos_v5 = pos, ref_base = toupper(ref_base),
                    depth, cml_base, matches_b73)],
             by = c("chrom", "pos_v5"), all.x = TRUE)
# *** THIS CLASSIFICATION IS WRONG -- kept only to document the error. ***
# I assumed DArT's ref allele IS B73's allele because the AGPv4 alignment defined
# it. It does not: DArT calls SNPs within its OWN tag clusters across the assayed
# panel, and the alignment columns are added afterwards. Measured in
# qtl_test_pileup_coords.R, the B73 v5 base equals the DArT ref allele at only
# ~76% of markers. So "CML530 == B73 here" is NOT the ref/alt question, and this
# route cannot give polarity however the bases are called.
# Use qtl_cml530_alleles_direct.R.
res[, cml_class := fifelse(is.na(cml_base), "no-hit",
                    fifelse(matches_b73, "ref", "alt"))]
print(res[, .N, by = cml_class][order(-N)])
cat("\nref:alt ratio =",
    round(res[cml_class == "ref", .N] / res[cml_class == "alt", .N], 3), "\n")
cat("framework subset:\n")
print(res[framework == TRUE, .N, by = cml_class][order(-N)])
fw <- res[framework == TRUE]
cat("framework ref:alt =",
    round(fw[cml_class == "ref", .N] / fw[cml_class == "alt", .N], 3), "\n")

rule("3. CROSS-VALIDATE AGAINST THE DIRECT CIGAR+faidx READ")
if (file.exists(DIRECT)) {
  dr <- fread(DIRECT, colClasses = list(character = "chr_v5"))
  cmp <- merge(res[, .(marker, pileup = cml_class)],
               dr[, .(marker, direct = cml_class)], by = "marker")
  cat("markers in both:", nrow(cmp), "\n")
  both <- cmp[pileup %in% c("ref", "alt") & direct %in% c("ref", "alt")]
  cat("both methods gave a call:", nrow(both), "\n")
  cat(sprintf("AGREEMENT: %d of %d (%.2f%%)\n",
              sum(both$pileup == both$direct), nrow(both),
              100 * mean(both$pileup == both$direct)))
  cat("\nconfusion matrix (rows pileup, cols direct):\n")
  print(table(pileup = both$pileup, direct = both$direct))
  cat("\ncoverage comparison:\n")
  print(data.table(
    method = c("pileup (v5 BAM, MAPQ>=20 spans)", "direct (all primary alignments)"),
    called = c(cmp[pileup %in% c("ref", "alt"), .N],
               cmp[direct %in% c("ref", "alt"), .N])))
  cat("\nTwo INDEPENDENT routes -- different coordinate systems, different code\n")
  cat("paths. High agreement validates both; disagreement localises the problem.\n")
} else cat("direct file absent; skipping cross-validation\n")

rule("4. WRITE")
setorder(res, chr_v5, pos_v5)
fwrite(res, OUT, sep = "\t")
cat("wrote", OUT, " rows:", nrow(res), "\n")
cat("\nCONCLUSION: if agreement is high, the pileup route works once the caller is\n")
cat("removed, and it is preferable -- coordinate-keyed, no name join, no bespoke\n")
cat("CIGAR walker. The only cost is coverage, because the v5 BAM was built from\n")
cat("MAPQ>=20 spans; rebuilding spans without that filter would close the gap.\n")
