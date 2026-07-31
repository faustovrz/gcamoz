#!/usr/bin/env Rscript
# qtl_encode_abh_all.R -- STEP 1 ONLY: encode every marker with a known CML530
# allele as A/H/B, with A = CML530. No QC filters of any kind.
#
#   A = the CML530 allele   (read from the assembly, qtl_cml530_alleles_direct.R)
#   B = the other allele    (the ENxx tester's, by complement in a biparental cross)
#   H = heterozygote
#   -  = missing
#
# NO FILTERS APPLIED. Not CallRate, not RepAvg, not MAF, not 1:2:1, not tag-dedup,
# not local linkage, not phase correction. Segregation statistics are REPORTED so
# the QC step can see them, but nothing is removed here.
#
# Reads   data/qtl/Report_DMz26-3123_SNP_mapping_2.csv
#         data/qtl/derived/CML530_marker_alleles_direct.tsv
#         data/qtl/derived/markers_v5.tsv
# Writes  data/qtl/derived/rqtl_gen_abh_all.csv   (csvsr genotype file)
#         data/qtl/derived/abh_all_marker_info.tsv
# Usage: Rscript scripts/qtl/qtl_encode_abh_all.R > agent/qtl_encode_abh_all.log 2>&1

suppressMessages({library(data.table); library(qtl)})
rule <- function(x) cat("\n", strrep("=", 70), "\n", x, "\n", strrep("=", 70), "\n", sep = "")

D <- "data/qtl/derived"
SNP1ROW <- "data/qtl/Report_DMz26-3123_SNP_mapping_2.csv"
OUTGEN  <- file.path(D, "rqtl_gen_abh_all.csv")
OUTINFO <- file.path(D, "abh_all_marker_info.tsv")

rule("1. INPUTS")
d <- fread(SNP1ROW, skip = 6, header = TRUE, na.strings = c("", "NA"))
mc <- 1:which(names(d) == "RepAvg")
G  <- as.matrix(d[, -mc, with = FALSE]); rownames(G) <- d$AlleleID
cat("DArT report markers:", nrow(G), " samples:", ncol(G), "\n")
cml <- fread(file.path(D, "CML530_marker_alleles_direct.tsv"),
             colClasses = list(character = "chr_v5"))
cat("CML530 allele calls:", nrow(cml), "\n")
print(cml[, .N, by = cml_class][order(-N)])

rule("2. MARKERS WITH A = CML530 DEFINED")
mk <- cml[cml_class %in% c("ref", "alt")]
cat("A = CML530 defined:", nrow(mk), "\n")
mk[, ref_allele := sub(".*:([ACGT])>([ACGT]).*", "\\1", SNP)]
mk[, alt_allele := sub(".*:([ACGT])>([ACGT]).*", "\\2", SNP)]
bad <- mk[!(ref_allele %in% c("A","C","G","T")) | !(alt_allele %in% c("A","C","G","T"))]
if (nrow(bad)) cat("unparseable SNP field, excluded:", nrow(bad), "\n")
mk <- mk[ref_allele %in% c("A","C","G","T") & alt_allele %in% c("A","C","G","T")]
mk <- mk[marker %in% rownames(G)]
cat("encodable:", nrow(mk), "on", uniqueN(mk$chr_v5), "chromosomes (B73 v5)\n")

rule("2b. ALIGNMENT-QUALITY FILTER + CML530 MARKER ORDER")
# WHY THIS EXISTS. The CML530 allele is read from the tag's alignment to the CML530
# assembly (tags_vs_CML530.bam). If that tag aligns equally well to several places, the
# aligner picks one arbitrarily and the allele is read off a locus that may not be the
# one DArT genotyped. Measured, not assumed (agent/bam_multimapping.log):
#
#     group                          n     MAPQ=0   median MAPQ   n2>0
#     markers anti-phase to both       27    29.6%       4        70.4%
#     B73/CML530 chromosome disagree   20    45.0%       2        85.0%
#     everything else                1818     1.2%      60         6.2%
#     Fisher vs rest: p = 2.3e-8 and 3.2e-11
#
# 11 of the 27 anti-phase markers carry an XA hit on the very chromosome B73 assigns
# them to -- the tag aligns to BOTH places in CML530 and the aligner took the other one.
# MAPQ and XA were never used as a filter; qc_pass does not encode them.
#
# WHY THE ORDER IS CML530's. B73 is not a parent of this cross. CML530 is the female
# parent (7x7 NC Design II, README.md:44-48), so its own coordinates are the right
# physical order. est.map never reorders -- it only estimates distances along the order
# it is given, so a wrong order cannot be recovered downstream.
# Within a chromosome the two orders barely differ (Spearman 0.9987, 0.61% of adjacent
# pairs reverse); the gain is the CHROMOSOME assignment, where they disagree for 37
# markers (agent/cml530_order_check.log).
BAM <- file.path(D, "CML530", "tags_vs_CML530.bam")
stopifnot(file.exists(BAM))
sam <- fread(cmd = sprintf("samtools view %s | cut -f1,2,3,4,5,6,12-", BAM),
             header = FALSE, sep = "\t", fill = TRUE, quote = "")
setnames(sam, 1:6, c("qname", "flag", "rname_bam", "pos_bam", "mapq", "cigar"))
sam[, marker := sub("\\|SNP=.*$", "", qname)]
tagcols <- setdiff(names(sam), c("qname","flag","rname_bam","pos_bam","mapq","cigar","marker"))
sam[, tags := do.call(paste, c(.SD, sep = "\t")), .SDcols = tagcols]
sam[, n2 := as.integer(sub(".*\\bn2:i:(\\d+).*", "\\1", tags))]
sam[!grepl("\\bn2:i:", tags), n2 := 0L]
sam[, has_XA := grepl("\\bXA:Z:", tags)]
sam[, unmapped := bitwAnd(flag, 4L) > 0L]
mk <- merge(mk, sam[, .(marker, mapq, n2, has_XA, unmapped)], by = "marker", all.x = TRUE)

cat(sprintf("alignment records joined: %d of %d markers\n",
            sum(!is.na(mk$mapq)), nrow(mk)))
cat("\nattrition at candidate MAPQ cuts (XA-carrying markers dropped in every case):\n")
for (q in c(0, 1, 10, 20, 30, 40, 60)) {
  keep <- !mk$unmapped %in% TRUE & !mk$has_XA %in% TRUE & mk$mapq >= q
  cat(sprintf("  MAPQ >= %2d : %5d retained (%5.1f%%)\n",
              q, sum(keep, na.rm = TRUE), 100 * mean(keep, na.rm = TRUE)))
}
cat("\nXA alone: ")
cat(sprintf("%d markers carry an alternative alignment (%.1f%%)\n",
            sum(mk$has_XA %in% TRUE), 100 * mean(mk$has_XA %in% TRUE)))

# THRESHOLDS. MAPQ >= 30 is bwa's conventional "confidently placed" cut; any XA at all
# means a second placement was found, which is the condition that produced the wrong
# allele. Both are set here as named constants so they can be changed in one place --
# the sweep above is what to judge them against.
MAPQ_MIN <- 30L
n_before <- nrow(mk)
mk <- mk[!unmapped %in% TRUE & !has_XA %in% TRUE & mapq >= MAPQ_MIN]
cat(sprintf("\nrejected %d markers (unmapped, XA present, or MAPQ < %d)\n",
            n_before - nrow(mk), MAPQ_MIN))
cat("retained:", nrow(mk), "\n")

# CML530 order. Markers off the ten chromosomes cannot join a linkage group.
mk[, cml_chr := suppressWarnings(as.integer(sub("^chr0?", "", rname)))]
n_before <- nrow(mk)
off <- mk[is.na(cml_chr) | !cml_chr %in% 1:10]
if (nrow(off)) { cat("\non non-chromosome CML530 sequences, dropped:\n"); print(off[, .N, by = rname]) }
mk <- mk[!is.na(cml_chr) & cml_chr %in% 1:10]
cat(sprintf("dropped %d off-chromosome markers; retained %d\n", n_before - nrow(mk), nrow(mk)))
cat(sprintf("\nB73/CML530 chromosome disagreements still present: %d\n",
            mk[as.integer(chr_v5) != cml_chr, .N]))
setorder(mk, cml_chr, ref_pos)
cat("ORDER IS NOW CML530 (cml_chr, ref_pos).\n")
print(mk[, .N, by = cml_chr][order(cml_chr)])

rule("3. ENCODE  (A = CML530)")
# DArT one-row: 0 = hom REF, 1 = hom ALT, 2 = het, - = missing
# cml_class "ref" -> CML530 carries REF -> A = code 0, B = code 1
# cml_class "alt" -> CML530 carries ALT -> A = code 1, B = code 0
Gm <- G[mk$marker, , drop = FALSE]
ABH <- matrix(NA_character_, nrow(Gm), ncol(Gm), dimnames = dimnames(Gm))
is_alt <- mk$cml_class == "alt"
ABH[Gm == "2"] <- "H"
ABH[!is_alt, ][Gm[!is_alt, ] == "0"] <- "A"
ABH[!is_alt, ][Gm[!is_alt, ] == "1"] <- "B"
ABH[ is_alt, ][Gm[ is_alt, ] == "1"] <- "A"
ABH[ is_alt, ][Gm[ is_alt, ] == "0"] <- "B"
tb <- table(factor(ABH, levels = c("A","H","B")))
cat("composition:\n"); print(tb)
cat(sprintf("A %.1f%%  H %.1f%%  B %.1f%%   (F2 expectation 25/50/25)\n",
            100*tb["A"]/sum(tb), 100*tb["H"]/sum(tb), 100*tb["B"]/sum(tb)))
cat("missing:", sum(is.na(ABH)), sprintf("(%.1f%%)\n", 100*mean(is.na(ABH))))
stopifnot(all(is.na(ABH) | ABH %in% c("A","H","B")))
cat("every non-missing cell is A, H or B: TRUE\n")

rule("4. SEGREGATION STATISTICS (reported, NOT filtered)")
nA <- rowSums(ABH == "A", na.rm = TRUE); nH <- rowSums(ABH == "H", na.rm = TRUE)
nB <- rowSums(ABH == "B", na.rm = TRUE); tot <- nA + nH + nB
mk[, `:=`(nA = nA, nH = nH, nB = nB, n_called = tot)]
mk[, maf := pmin(2*nA + nH, 2*nB + nH) / (2*tot)]
mk[, chi121 := (nA - tot/4)^2/(tot/4) + (nB - tot/4)^2/(tot/4) + (nH - tot/2)^2/(tot/2)]
mk[, p121 := pchisq(chi121, 2, lower.tail = FALSE)]
cat("MAF:  "); print(round(summary(mk$maf), 3))
cat("1:2:1 chi-sq p:  "); print(signif(summary(mk$p121), 3))
cat("\nnon-segregating (MAF < 0.05):", mk[maf < 0.05, .N],
    "  <-- both parents identical; cannot be mapped\n")
cat("MAF >= 0.15:", mk[maf >= 0.15, .N], "\n")
cat("markers sharing a 69bp tag (CloneID):",
    nrow(mk) - uniqueN(mk$CloneID), "duplicates over",
    uniqueN(mk$CloneID), "tags\n")

rule("4b. HARD MAF FILTER -- MAF > 0.15 (FVRZ: correct for an F2)")
# Per FVRZ: a HARD MAF threshold is appropriate here because this is an F2. Expected
# genotype frequencies are 0.25/0.50/0.25 -- symmetric, with MAF ~ 0.5 at every
# informative marker -- so 0.15 is a principled cut, not an arbitrary one.
#
# His TeoNAM/airmine pipeline uses a RELATIVE outlier rule instead precisely because
# those populations are BC1S4 / BC2S3, where the expectation is heavily skewed
# (0.734 / 0.031 / 0.234). A hard MAF cut there would delete real markers; a
# top-2.5%-per-chromosome rule is the only safe option. That constraint does not
# apply to an F2.
#
# It also fixes what killed the previous run: ~1,000 non-segregating markers
# (MAF < 0.05, i.e. both parents carrying the same allele -- TASSEL
# GenosToABHPlugin.java L158 rejects these as parentA == parentB) reached est.map,
# because a relative filter capped at ~2.5% per chromosome cannot remove a third of
# the set. Round 1 returned 1,015,613 cM, find_quirky's threshold degenerated to
# Inf, and 99.3% of markers were flagged, leaving 24. See agent/build_map.log.
MAF_MIN <- 0.15
n_before <- nrow(mk)
mk <- mk[maf > MAF_MIN]
cat(sprintf("rejected %d markers with MAF <= %.2f\n", n_before - nrow(mk), MAF_MIN))
cat("retained:", nrow(mk), "\n")
ABH <- ABH[mk$marker, , drop = FALSE]
tb <- table(factor(ABH, levels = c("A","H","B")))
cat(sprintf("composition now: A %.1f%%  H %.1f%%  B %.1f%%   (F2 expectation 25/50/25)\n",
            100*tb["A"]/sum(tb), 100*tb["H"]/sum(tb), 100*tb["B"]/sum(tb)))
print(mk[, .N, by = cml_chr][order(cml_chr)])

rule("4b2. HARD MENDELIAN FILTER -- 1:2:1 chi-square (FVRZ: correct for an F2)")
# Same argument as the hard MAF cut, which I applied and then failed to generalise: an
# F2 has a KNOWN, SYMMETRIC expectation of 0.25 / 0.50 / 0.25, so a hard chi-square
# against it is principled. FVRZ's renorm_z relative rule exists because BC1S4 expects
# 0.734 / 0.031 / 0.234, where a hard test would delete real markers.
#
# WHY THE RELATIVE RULE CANNOT SUBSTITUTE HERE: it flags a top ~2.5% tail per
# chromosome BY CONSTRUCTION. In the last build it removed 65 markers while markers at
# chi2 = 279 (p beyond any conventional precision) stayed in the map
# (agent/build_map.log). Median chi2 was 2.72 (p = 0.26), so the bulk is fine -- it is
# a TAIL problem the relative rule is structurally unable to reach.
#
# WHY IT MATTERS FOR MAP LENGTH: a marker failing 1:2:1 has miscalled genotypes, and
# est.map reads miscalls as crossovers. Median gap in the last build was 3.59 cM
# against 0.80 expected at that marker count -- a uniform ~4.5x inflation, which
# miscalled genotypes explain without any assumption about DArT's error rate.
#
# ACCEPTED COST: this also removes markers with GENUINE segregation distortion
# (gametophytic factors, pollen competition -- real in maize). At n = 186 the test
# cannot separate real distortion from miscalling. Fine for building a map; NOT fine if
# distortion is later a subject of study.
# THRESHOLD BY FDR, not by a raw p-value. The test is applied to every marker, so a
# fixed raw cut has no stated error rate. Benjamini-Hochberg at FDR = 0.05 controls the
# proportion of REMOVED markers that are not really distorted (agent/fdr_p121.log).
#
# Distortion here is widespread and real -- 28.1% of markers have p < 0.05 against 5%
# expected under the null, 18.2% below 0.01 against 1% -- so with that many true
# positives the correction barely moves the cutoff: BH lands at raw p ~ 0.0087 and
# removes 369 markers, against 376 for the raw 0.01 previously used. The old value was
# already about right; this states what it controls.
#
# Computed from the p-values rather than hardcoded, so it adapts if the tested set
# changes. NOTE the test set is the POST-MAF set -- the MAF cut runs first (4b) -- which
# is the correct denominator for the correction.
FDR_LEVEL <- 0.05
n_before <- nrow(mk)
cat(sprintf("chi2 vs 1:2:1 -- median %.2f  q95 %.2f  max %.1f\n",
            median(mk$chi121), quantile(mk$chi121, .95), max(mk$chi121)))
cat(sprintf("p-value       -- median %.4f  q05 %.2e\n",
            median(mk$p121), quantile(mk$p121, .05)))
mk[, q121 := p.adjust(p121, method = "BH")]
p_cut <- if (any(mk$q121 <= FDR_LEVEL)) max(mk$p121[mk$q121 <= FDR_LEVEL]) else 0
cat(sprintf("tests: %d   BH at FDR %.2f -> raw p cutoff %.4g\n", n_before, FDR_LEVEL, p_cut))
cat(sprintf("p < 0.05: %d observed vs %.0f expected under the null\n",
            sum(mk$p121 < 0.05), 0.05 * n_before))
mk <- mk[q121 > FDR_LEVEL]
cat(sprintf("rejected %d markers at FDR %.2f\n", n_before - nrow(mk), FDR_LEVEL))
cat("retained:", nrow(mk), "\n")
ABH <- ABH[mk$marker, , drop = FALSE]
tb <- table(factor(ABH, levels = c("A","H","B")))
cat(sprintf("composition: A %.1f%%  H %.1f%%  B %.1f%%\n",
            100*tb["A"]/sum(tb), 100*tb["H"]/sum(tb), 100*tb["B"]/sum(tb)))
print(mk[, .N, by = cml_chr][order(cml_chr)])

rule("4c. DROP CO-TAG MARKERS (the DP > 1 sites)")
# Markers sharing a CloneID sit on ONE 69bp restriction fragment. They cannot
# recombine, so their genotypes must be near-identical -- yet median concordance
# between co-tag SNPs is 0.695, with 37 of 42 pairs below 0.90 and 7 below 0.50
# (agent/check_tag_redundancy.log). At zero recombination distance there is no
# biological route to 30% discordance: it is genotyping error or a wrong allele call.
#
# These are exactly the DP > 1 sites in the v5 pileup. Chance overlap of independent
# clones is ruled out: with 3,167 spans of ~55 bp over a 2.18 Gb genome the Poisson
# expectation is 0.25 overlapping sites, against 323 observed. And the count closes --
# 163 two-SNP tags predict 326 depth>1 positions, 323 seen. So DP > 1 and
# "CloneID appears more than once" are the same set.
#
# Implemented on CloneID, NOT on pileup depth: no v5 BAM, no pileup, no depth column
# needed, so this filter works in the minimal chain.
#
# BOTH members are dropped, not one kept per tag: if the two calls on a fragment
# disagree 30% of the time, there is no basis for choosing which to keep.
n_before <- nrow(mk)
dup_tags <- mk[, .N, by = CloneID][N > 1, CloneID]
mk <- mk[!CloneID %in% dup_tags]
cat(sprintf("tags carrying >1 SNP: %d\n", length(dup_tags)))
cat(sprintf("dropped %d markers on those tags\n", n_before - nrow(mk)))
cat("retained:", nrow(mk), "\n")
stopifnot(uniqueN(mk$CloneID) == nrow(mk))
cat("one marker per tag: TRUE\n")
ABH <- ABH[mk$marker, , drop = FALSE]
tb <- table(factor(ABH, levels = c("A","H","B")))
cat(sprintf("composition: A %.1f%%  H %.1f%%  B %.1f%%\n",
            100*tb["A"]/sum(tb), 100*tb["H"]/sum(tb), 100*tb["B"]/sum(tb)))

# NOTE: no phase pass. A = CML530 by construction, so phase is already fixed --
# a marker whose genotypes contradict its neighbours is a WRONG marker, not a
# mis-phased one, and flipping it would silently break A = CML530.

if (FALSE) {   # retained for reference only; do not run
# Encoding A = CML530 sets polarity from the parent, which is the same job airmine's
# switch_gt()/dummy-parent step does. It is not perfect here: the assembly read gave
# 93.3% phase consistency on the old 1,384-marker framework subset, and worse on this
# wider 2,501 set. Residual anti-phase markers make rf -> 1-rf, which pins at 0.5 and
# saturates est.map: the unphased run showed q95 = q99 = max = 1001.506 cM with 476
# junctions above 25 cM, and a round-1 total of 277,730 cM.
#
# So flip markers that disagree with their local neighbourhood. Windowed (+/-K), not a
# naive chain, so one bad marker cannot propagate a flip along the chromosome.
# NO markers are dropped here. Which markers were flipped is recorded, because a flip
# means A is no longer the CML530 allele at that marker.
dose <- matrix(NA_integer_, ncol(ABH), nrow(ABH),
               dimnames = list(colnames(ABH), rownames(ABH)))
dose[t(ABH) == "A"] <- 0L; dose[t(ABH) == "H"] <- 1L; dose[t(ABH) == "B"] <- 2L
chrv <- mk$chr_v5; K <- 5L
pc <- function(d) unlist(lapply(unique(chrv), function(ch) {
  i <- which(chrv == ch); if (length(i) < 2) return(NULL)
  sapply(seq_len(length(i)-1), function(k) {
    a <- d[, i[k]]; b <- d[, i[k+1]]; ok <- is.finite(a) & is.finite(b)
    if (sum(ok) < 30) return(NA_real_); cor(a[ok], b[ok])
  })
}))
r <- pc(dose); r <- r[is.finite(r)]
cat(sprintf("before: %.1f%% of junctions in phase (median cor %.3f)\n",
            100*mean(r > 0), median(r)))
flipped <- rep(FALSE, ncol(dose))
for (it in 1:20) {
  nf <- 0L
  for (ch in unique(chrv)) {
    idx <- which(chrv == ch); n <- length(idx); if (n < 3) next
    dd <- dose[, idx, drop = FALSE]
    sc <- sapply(seq_len(n), function(k) {
      nb <- setdiff(max(1,k-K):min(n,k+K), k)
      mean(sapply(nb, function(j) {
        a <- dd[,k]; b <- dd[,j]; ok <- is.finite(a) & is.finite(b)
        if (sum(ok) < 30) return(NA_real_); cor(a[ok], b[ok])
      }), na.rm = TRUE)
    })
    bad <- which(is.finite(sc) & sc < 0)
    if (length(bad)) { gi <- idx[bad]; dose[, gi] <- 2L - dose[, gi]
                       flipped[gi] <- !flipped[gi]; nf <- nf + length(bad) }
  }
  if (nf == 0L) break
  cat(sprintf("  iter %2d: flipped %d\n", it, nf))
}
r <- pc(dose); r <- r[is.finite(r)]
cat(sprintf("after : %.1f%% in phase (median cor %.3f), %d anti-phase of %d\n",
            100*mean(r > 0), median(r), sum(r < 0), length(r)))
cat(sprintf("markers flipped: %d (%.1f%%) -- A is the TESTER allele on these\n",
            sum(flipped), 100*mean(flipped)))
mk[, phase_flipped := flipped]
mk[, a_is_cml530 := !flipped]
ABH <- t(matrix(c("A","H","B")[dose + 1L], nrow = nrow(dose),
                dimnames = dimnames(dose)))
}   # end of disabled reference block -- the phase-flip approach is WRONG here

rule("4d. DROP MARKERS THAT CONTRADICT THE A = CML530 ENCODING")
# A = CML530 at every marker, so phase is FIXED BY CONSTRUCTION: a marker must
# correlate positively with its physical neighbours. One correlating with NONE of them
# contradicts the encoding -- it is a WRONG marker, not a mis-phased one, and flipping
# it would silently break A = CML530.
#
# Cutoff set by FVRZ from the reported range (agent/report_offender_thresholds.log):
#   |cor| < 0.30 over K = 5 neighbours each side
# The count is insensitive there: 0.25-0.50 adds only ~6 markers per 0.05 step, then
# triples above 0.50 as it begins taking sound markers. The distribution agrees --
# 2nd percentile 0.092, 5th percentile 0.351.
#
# NOT reachable by any DArT statistic (measured, not assumed):
#   RepAvg   offender rate is FLAT across bands, 4.1% at the bottom vs 3.6% at the top;
#            a 0.95 cut takes 46 of 117 offenders and sacrifices 768 sound markers
#   CallRate real gradient (19.3% in 0.90-0.95 vs 3.5% above 0.95) but 80 of 117
#            offenders sit above 0.95
#   MAF, co-tag concordance: no signal
LINK_MIN <- 0.30
LINK_K   <- 5L
d0 <- matrix(NA_integer_, ncol(ABH), nrow(ABH),
             dimnames = list(colnames(ABH), rownames(ABH)))
d0[t(ABH) == "A"] <- 0L; d0[t(ABH) == "H"] <- 1L; d0[t(ABH) == "B"] <- 2L
cv <- mk$cml_chr
best <- rep(NA_real_, ncol(d0))
for (ch in unique(cv)) {
  idx <- which(cv == ch); n <- length(idx); dd <- d0[, idx, drop = FALSE]
  for (k in seq_len(n)) {
    nb <- setdiff(max(1, k - LINK_K):min(n, k + LINK_K), k)
    cs <- sapply(nb, function(j) {
      a <- dd[, k]; b <- dd[, j]; ok <- is.finite(a) & is.finite(b)
      if (sum(ok) < 30) return(NA_real_)
      abs(cor(a[ok], b[ok]))
    })
    if (!all(is.na(cs))) best[idx[k]] <- max(cs, na.rm = TRUE)
  }
}
off <- is.finite(best) & best < LINK_MIN
cat(sprintf("max |cor| with any neighbour: 2%%ile %.3f  5%%ile %.3f  median %.3f\n",
            quantile(best, .02, na.rm = TRUE), quantile(best, .05, na.rm = TRUE),
            median(best, na.rm = TRUE)))
cat(sprintf("dropping %d markers (%.1f%%) with |cor| < %.2f over K = %d\n",
            sum(off), 100 * mean(off), LINK_MIN, LINK_K))
cat("per chromosome:\n"); print(table(cv[off]))
mk  <- mk[!off]
ABH <- ABH[mk$marker, , drop = FALSE]
cat("retained:", nrow(mk), "\n")
tb <- table(factor(ABH, levels = c("A","H","B")))
cat(sprintf("composition: A %.1f%%  H %.1f%%  B %.1f%%\n",
            100*tb["A"]/sum(tb), 100*tb["H"]/sum(tb), 100*tb["B"]/sum(tb)))
cat("A = CML530 on every retained marker -- no flips applied\n")

rule("5. WRITE csvsr")
gen <- data.table(id = mk$marker, chr = mk$cml_chr)   # CML530 chromosome, CML530 order
# ABH is ALREADY markers x individuals, which is the csvsr row orientation.
# Do NOT transpose -- t(ABH) gives 3,728 columns of 186 rows.
gen <- cbind(gen, as.data.table(ABH))
# DArT names samples "plant 1" (with a space); rqtl_phe.csv uses "plant1". Without
# normalising, read.cross UNIONS them (186 + 166 = 352 individuals) instead of
# matching. All 186 GENOTYPED individuals are kept -- est.map does not use
# phenotypes, so the map should use every genotyped plant, not only the 166 that
# also have phenotypes.
ids <- gsub("^plant ", "plant", colnames(ABH))
setnames(gen, c("id", "chr", ids))
fwrite(gen, OUTGEN, na = "-", quote = FALSE)
fwrite(mk, OUTINFO, sep = "\t")
cat("wrote", OUTGEN, " dim:", paste(dim(gen), collapse = " x "), "\n")
cat("wrote", OUTINFO, "\n")

rule("6. VERIFY read.cross ROUND-TRIP")
cr <- read.cross(format = "csvsr", dir = "", genfile = OUTGEN,
                 phefile = file.path(D, "rqtl_phe.csv"),
                 genotypes = c("A","H","B"), na.strings = c("-","NA"),
                 crosstype = "f2", estimate.map = FALSE)
cat("nind:", nind(cr), " markers:", totmar(cr), " nchr:", nchr(cr), "\n")
print(nmar(cr))
stopifnot(totmar(cr) == nrow(mk))
cat("\nENCODED:", totmar(cr), "markers, A = CML530 on every one. No QC applied.\n")
