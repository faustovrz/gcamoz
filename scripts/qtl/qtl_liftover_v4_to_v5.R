#!/usr/bin/env Rscript
# qtl_liftover_v4_to_v5.R -- lift the DArTseq SNP positions for order DMz26-3123
# from B73 AGPv4 (the assembly DArT aligned against) to B73 NAM v5.
#
# WHY: MaizeGDB annotation and the candidate-gene literature are on v5, so Q5
# (candidate genes under QTL peaks) needs v5 coordinates. Kept deliberately
# MAP-NEUTRAL -- this produces a marker roster + positions only, no cM, matching
# the convention of zealhmm's data/teonam/markers_v5.tsv.
#
# Reuses lift_unique() from ~/repos/zealhmm/R/teonam_liftover.R verbatim so this
# behaves identically to the TeoNAM v2->v4->v5 and molbreeding v3->v4->v5 lifts:
# ONE hop, keep only unique 1:1 mappings.
#
# Lifts ALL v4-anchored markers, not just the QC-passing set, so the liftover is
# independent of any marker-filtering decision; framework membership is recorded
# as a column instead.
#
# Reads   data/qtl/Report_DMz26-3123_SNP_mapping_2.csv   (read-only)
# Writes  data/qtl/derived/markers_v5.tsv                (roster + v4/v5 coords)
#         data/qtl/derived/markers_v5_qc.csv             (incl. failures + reasons)
#         data/qtl/derived/tags_for_cml530.fa            (tags, for allele calling)
# stdout -> agent/qtl_liftover_v4_to_v5.log
#
# Run: Rscript agent/qtl_liftover_v4_to_v5.R > agent/qtl_liftover_v4_to_v5.log 2>&1

suppressMessages({
  library(data.table); library(GenomicRanges); library(IRanges)
  library(rtracklayer); library(S4Vectors); library(Biostrings)
})

SNP1ROW   <- "data/qtl/Report_DMz26-3123_SNP_mapping_2.csv"
DART_SKIP <- 6
CHAIN     <- "/Users/fvrodriguez/repos/zealhmm/data/ref/chain_files/B73_RefGen_v4_to_Zm-B73-REFERENCE-NAM-5.0.chain"
OUTDIR    <- "data/qtl/derived"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
rule <- function(x) cat("\n", strrep("=", 70), "\n", x, "\n", strrep("=", 70), "\n", sep = "")

# ---- house liftover helper (verbatim from zealhmm/R/teonam_liftover.R) --------
lift_unique <- function(gr, chain_path) {
  stopifnot(file.exists(chain_path))
  ch <- rtracklayer::import.chain(chain_path)
  lifted <- rtracklayer::liftOver(gr, ch) # GRangesList parallel to gr
  n1 <- S4Vectors::elementNROWS(lifted) == 1L # exactly one target range
  out <- unlist(lifted[n1])
  out$marker <- gr$marker[n1]
  out
}

## ------------------------------------------------------------------ 1. input
rule("1. INPUT")
d <- fread(SNP1ROW, skip = DART_SKIP, header = TRUE, na.strings = c("", "NA"))
mc <- 1:which(names(d) == "RepAvg")
G  <- as.matrix(d[, -mc, with = FALSE])
setnames(d,
  c("Chrom_Maize_B73_V4.0.assembly", "ChromPosSnp_Maize_B73_V4.0.assembly",
    "ChromPosTag_Maize_B73_V4.0.assembly", "AlnCnt_Maize_B73_V4.0.assembly",
    "Strand_Maize_B73_V4.0.assembly"),
  c("chrom_raw", "pos_snp_v4", "pos_tag_v4", "aln_cnt", "strand_v4"))

d[, chr_v4 := sub(" .*", "", chrom_raw)]
d[, pos_v4 := suppressWarnings(as.numeric(pos_snp_v4))]
cat("markers in report:", nrow(d), " samples:", ncol(G), "\n")
cat("v4-anchored to chr 1-10:", sum(d$chr_v4 %in% as.character(1:10)), "\n")
cat("chain:", basename(CHAIN), "\n")

# recompute the QC/framework flags (same rule as agent/qtl_data_audit.R) so the
# roster can record framework membership without re-deriving it downstream
n0 <- rowSums(G == "0"); n1 <- rowSums(G == "1"); n2 <- rowSums(G == "2")
tot <- n0 + n1 + n2
pref <- (2 * n0 + n2) / (2 * tot); maf <- pmin(pref, 1 - pref)
cs <- (n0 - tot/4)^2/(tot/4) + (n1 - tot/4)^2/(tot/4) + (n2 - tot/2)^2/(tot/2)
p121 <- pchisq(cs, 2, lower.tail = FALSE)
d[, `:=`(maf = maf, p121 = p121)]
d[, qc_pass := chr_v4 %in% as.character(1:10) & CallRate >= 0.90 &
     RepAvg >= 0.95 & maf >= 0.15 & p121 > 0.01]
d[is.na(qc_pass), qc_pass := FALSE]
# framework = 1 SNP per CloneID (highest CallRate) among qc_pass
d[, framework := FALSE]
fwi <- d[qc_pass == TRUE][order(-CallRate), .SD[1], by = CloneID]$AlleleID
d[AlleleID %in% fwi, framework := TRUE]
cat("qc_pass:", sum(d$qc_pass), " framework (1 SNP/tag):", sum(d$framework), "\n")

## --------------------------------------------------------------- 2. liftover
rule("2. LIFTOVER v4 -> v5 (single hop, unique 1:1 only)")
ok <- d$chr_v4 %in% as.character(1:10) & is.finite(d$pos_v4)
gr <- GRanges(seqnames = d$chr_v4[ok],
              ranges   = IRanges(start = as.integer(d$pos_v4[ok]), width = 1L),
              marker   = d$AlleleID[ok])
cat("submitted to liftOver:", length(gr), "\n")
gr5 <- lift_unique(gr, CHAIN)
cat("lifted 1:1            :", length(gr5),
    sprintf(" (%.1f%%)\n", 100 * length(gr5) / length(gr)))
cat("failed / multi-mapping:", length(gr) - length(gr5), "\n")

lift <- data.table(marker = gr5$marker,
                   chr_v5 = as.character(seqnames(gr5)),
                   pos_v5 = start(gr5))
res <- merge(d[, .(marker = AlleleID, CloneID, chr_v4, pos_v4, SNP, SnpPosition,
                   strand_v4, aln_cnt, CallRate, RepAvg, maf, p121,
                   qc_pass, framework,
                   TrimmedSequenceRef, TrimmedSequenceSnp)],
             lift, by = "marker", all.x = TRUE)
res[, lifted := !is.na(pos_v5)]

## ------------------------------------------------------- 3. QC on the lift
rule("3. LIFTOVER QC")
res[, off_chr10 := lifted & !(chr_v5 %in% as.character(1:10))]
res[, chr_change := lifted & chr_v5 != chr_v4]
cat("lifted onto a non-chr1-10 scaffold:", sum(res$off_chr10, na.rm = TRUE), "\n")
cat("CHROMOSOME CHANGERS (v4 chr != v5 chr):", sum(res$chr_change, na.rm = TRUE), "\n")
if (sum(res$chr_change, na.rm = TRUE)) {
  cat("  by v4 chr -> v5 chr:\n")
  print(res[chr_change == TRUE, .N, by = .(chr_v4, chr_v5)][order(-N)][1:min(15, .N)])
}
cat("\nNOTE: chromosome-changers are dropped from the roster (house rule, cf.\n",
    "liftover_teonam(): res$chr_v5 == res$chr_v2). They are the same class of\n",
    "artefact that the chr7 v2->v5 displaced block was, and are exactly what\n",
    "find_quirky()'s island rule is there to catch if any slip through.\n")

# keep = lifted, on chr 1-10, same chromosome
res[, keep := lifted & !off_chr10 & !chr_change]
res[is.na(keep), keep := FALSE]
cat("\nroster (lifted, chr1-10, same chr):", sum(res$keep), "\n")

# order inversions: rank disagreement between v4 and v5 within a chromosome
r <- res[keep == TRUE][order(as.integer(chr_v5), pos_v5)]
r[, rank_v5 := seq_len(.N), by = chr_v5]
r[, rank_v4 := frank(pos_v4, ties.method = "first"), by = chr_v4]
r[, inversion := rank_v5 != rank_v4]
cat("local order changes v4->v5 (rank mismatch):", sum(r$inversion), "\n")
cat("per-chromosome Spearman rho (v4 pos vs v5 pos):\n")
print(r[, .(n = .N, rho = round(cor(pos_v4, pos_v5, method = "spearman"), 5)),
        by = chr_v5][order(as.integer(chr_v5))])

cat("\nretention by class:\n")
print(res[, .(n = .N, lifted = sum(lifted), kept = sum(keep),
              pct_kept = round(100 * sum(keep) / .N, 1)),
          by = .(anchored_v4 = chr_v4 %in% as.character(1:10), qc_pass, framework)][order(-n)])

## ------------------------------------------------------------- 4. framework
rule("4. FRAMEWORK MAP IN v5 COORDINATES")
fw <- r[framework == TRUE]
cat("framework markers surviving the lift:", nrow(fw), "of", sum(d$framework), "\n")
fw[, mb := pos_v5 / 1e6]
print(fw[order(as.integer(chr_v5), pos_v5),
         .(n = .N, span_Mb = round(max(mb) - min(mb), 1),
           median_gap_Mb = round(median(diff(mb)), 2),
           max_gap_Mb = round(max(diff(mb)), 1)),
         by = chr_v5][order(as.integer(chr_v5))])

## --------------------------------------------------------------- 5. outputs
rule("5. OUTPUTS")
roster <- r[, .(marker, CloneID, chr_v4, pos_v4, chr_v5, pos_v5, rank_v5,
                inversion, SNP, SnpPosition, strand_v4,
                CallRate, RepAvg, maf, p121, qc_pass, framework)]
setorder(roster, chr_v5, pos_v5)
fwrite(roster, file.path(OUTDIR, "markers_v5.tsv"), sep = "\t")
fwrite(res[, .(marker, CloneID, chr_v4, pos_v4, chr_v5, pos_v5,
               lifted, off_chr10, chr_change, keep, qc_pass, framework)],
       file.path(OUTDIR, "markers_v5_qc.csv"))

# tag FASTA for the CML530 allele-calling step. TrimmedSequenceRef carries the
# REFERENCE allele of the SNP; SnpPosition is the 0-based offset within the tag.
tags <- res[keep == TRUE & !is.na(TrimmedSequenceRef) & TrimmedSequenceRef != ""]
fa <- DNAStringSet(tags$TrimmedSequenceRef)
names(fa) <- sprintf("%s|SNP=%s|snppos=%s|v5=%s:%s",
                     tags$marker, tags$SNP, tags$SnpPosition,
                     tags$chr_v5, tags$pos_v5)
writeXStringSet(fa, file.path(OUTDIR, "tags_for_cml530.fa"))
cat("wrote", file.path(OUTDIR, "markers_v5.tsv"), "  rows:", nrow(roster), "\n")
cat("wrote", file.path(OUTDIR, "markers_v5_qc.csv"), " rows:", nrow(res), "\n")
cat("wrote", file.path(OUTDIR, "tags_for_cml530.fa"), " seqs:", length(fa),
    " width:", paste(range(width(fa)), collapse = "-"), "bp\n")
cat("\nliftover complete.\n")
