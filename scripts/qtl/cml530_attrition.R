#!/usr/bin/env Rscript
# qtl_cml530_attrition.R -- reconcile the CML530 pipeline funnel EXACTLY.
#
# WHY: I reported "3740 mapped (97.7%)" and "651 MAPQ<20" and "no-hit is 708" as
# if they were independent, which does not add up. `samtools view -F 0x904`
# strips flag 0x4 (UNMAPPED) as well as 0x100 (secondary) and 0x800
# (supplementary), so the 651 records removed before span extraction ALREADY
# CONTAIN the 88 unmapped tags. Reporting 88 and 651 as separate losses
# double-counts the 88.
#
# This script derives every stage count from the artefacts and asserts that the
# losses sum to the observed no-hit total, so the funnel is arithmetic rather
# than narration.
#
# READ-ONLY.
# Usage: Rscript agent/qtl_cml530_attrition.R > agent/qtl_cml530_attrition.log 2>&1

suppressMessages(library(data.table))
rule <- function(x) cat("\n", strrep("=", 70), "\n", x, "\n", strrep("=", 70), "\n", sep = "")

D    <- "data/qtl/derived"
OUT  <- file.path(D, "CML530")
BAM1 <- file.path(OUT, "tags_vs_CML530.bam")
BAM2 <- file.path(OUT, "CML530_vs_B73v5.bam")
sam  <- function(...) as.integer(system2("samtools", c("view", "-c", ...), stdout = TRUE))

rule("1. STAGE COUNTS (from the BAMs, not from memory)")
n_tags      <- sam(BAM1)                         # every tag has a record
n_unmapped  <- sam("-f", "4", BAM1)
n_mapped    <- sam("-F", "4", BAM1)
n_secsup    <- sam("-f", "2304", BAM1)           # 0x100 | 0x800
n_kept      <- sam("-F", "0x904", "-q", "20", BAM1)
n_lowq      <- sam("-F", "0x904", BAM1) - n_kept # mapped, primary, MAPQ<20
n_span      <- length(readLines(file.path(OUT, "CML530_spans.regions")))
n_v5_in     <- sam(BAM2)
n_v5_mapped <- sam("-F", "4", BAM2)
n_v5_unmap  <- sam("-f", "4", BAM2)
n_vcf       <- length(system2("bcftools", c("view", "-H",
                 file.path(OUT, "CML530_at_markers.vcf")), stdout = TRUE))

cat(sprintf("tags in BAM1                        : %5d\n", n_tags))
cat(sprintf("  unmapped (flag 0x4)               : %5d\n", n_unmapped))
cat(sprintf("  mapped                            : %5d  (%.2f%%)\n", n_mapped,
            100 * n_mapped / n_tags))
cat(sprintf("  secondary/supplementary           : %5d\n", n_secsup))
cat(sprintf("  mapped primary but MAPQ<20        : %5d\n", n_lowq))
cat(sprintf("spans extracted (-F 0x904 -q 20)    : %5d\n", n_kept))
cat(sprintf("  regions file lines                : %5d\n", n_span))
cat(sprintf("spans submitted to B73 v5           : %5d\n", n_v5_in))
cat(sprintf("  mapped to v5                      : %5d  (%.2f%%)\n", n_v5_mapped,
            100 * n_v5_mapped / n_v5_in))
cat(sprintf("  failed to map to v5               : %5d\n", n_v5_unmap))
cat(sprintf("VCF records at marker positions     : %5d\n", n_vcf))

rule("2. WHAT '651' ACTUALLY WAS")
cat(sprintf("n_tags - n_spans = %d - %d = %d\n", n_tags, n_kept, n_tags - n_kept))
cat(sprintf("  = unmapped (%d) + secondary/supplementary (%d) + MAPQ<20 (%d) = %d\n",
            n_unmapped, n_secsup, n_lowq, n_unmapped + n_secsup + n_lowq))
stopifnot(n_tags - n_kept == n_unmapped + n_secsup + n_lowq)
cat("=> the 651 CONTAINS the 88 unmapped. They are not additive.\n")

rule("3. RECONCILE AGAINST no-hit")
al <- fread(file.path(D, "CML530_marker_alleles.tsv"),
            colClasses = list(character = "chr_v5"))
obs <- al[, .N, by = cml_class][order(-N)]
print(obs)
n_nohit  <- al[cml_class == "no-hit", .N]
n_called <- nrow(al) - n_nohit

# markers sharing one v5 position: a single VCF record serves >1 marker, so
# markers-with-a-call EXCEEDS the VCF record count by the number of extra sharers
dupes <- al[, .N, by = .(chr_v5, pos_v5)][N > 1]
n_extra <- sum(dupes$N) - nrow(dupes)
cat(sprintf("\nroster markers            : %5d\n", nrow(al)))
cat(sprintf("markers with a call       : %5d\n", n_called))
cat(sprintf("VCF records               : %5d\n", n_vcf))
cat(sprintf("markers sharing a v5 pos  : %5d extra sharers over %d positions\n",
            n_extra, nrow(dupes)))
cat(sprintf("  %d VCF records + %d extra sharers = %d markers with a call\n",
            n_vcf, n_extra, n_vcf + n_extra))
stopifnot(n_called == n_vcf + n_extra)

rule("4. THE FUNNEL, LOSSES SUMMED")
loss <- data.table(
  stage = c("no CML530 locus (unmapped)",
            "secondary/supplementary",
            "mapped to CML530 but MAPQ<20",
            "span failed to map back to B73 v5",
            "mapped to v5 but no VCF record at the marker position"),
  n = c(n_unmapped, n_secsup, n_lowq, n_v5_unmap, n_v5_mapped - n_vcf))
loss <- rbind(loss, data.table(stage = "TOTAL LOSSES", n = sum(loss$n)))
loss <- rbind(loss, data.table(
  stage = "less markers sharing a v5 position (recovered)", n = -n_extra))
loss <- rbind(loss, data.table(stage = "= no-hit", n = sum(loss$n[c(1:5, 7)])))
print(loss)
cat(sprintf("\nobserved no-hit in the classification: %d\n", n_nohit))
stopifnot(sum(loss$n[c(1:5, 7)]) == n_nohit)
cat("RECONCILED.\n")

rule("5. CORRECTED STATEMENT")
cat(sprintf("Of %d tags: %d (%.1f%%) found a CML530 locus; %d were dropped for\n",
            n_tags, n_mapped, 100 * n_mapped / n_tags, n_lowq))
cat(sprintf("MAPQ<20; %d spans failed to place back on B73 v5; %d placed but\n",
            n_v5_unmap, n_v5_mapped - n_vcf))
cat(sprintf("yielded no pileup record. Net: %d markers unpolarized (no-hit),\n", n_nohit))
cat(sprintf("%d polarized from the actual CML530 allele.\n", n_called))
cat(sprintf("\nFramework subset: %d of %d polarized (%.1f%%).\n",
            al[framework == TRUE & cml_class != "no-hit", .N],
            al[framework == TRUE, .N],
            100 * al[framework == TRUE & cml_class != "no-hit", .N] /
                  al[framework == TRUE, .N]))
cat("\nThe single largest loss is the MAPQ<20 filter, NOT absence of a CML530\n")
cat("locus. Those markers are still fully usable for MAPPING (recombination is\n")
cat("orientation-invariant); only their allele-effect SIGN is unavailable.\n")
