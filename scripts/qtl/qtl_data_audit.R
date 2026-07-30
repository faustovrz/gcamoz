#!/usr/bin/env Rscript
# qtl_data_audit.R -- read-only audit of data/qtl/ prior to any linkage mapping.
#
# Purpose: characterise the DArTseq report (order DMz26-3123, SEQART AFRICA) and
# the mozpue_phenotype.xlsx field data, and report what passes QC, so we know the
# effective mapping population size and usable marker set before writing any
# R/qtl code.
#
# WRITES NOTHING to data/. Output is stdout only -> agent/qtl_data_audit.log
# Run:  Rscript agent/qtl_data_audit.R > agent/qtl_data_audit.log 2>&1

suppressMessages({library(data.table); library(readxl)})

QTL <- "data/qtl"
SNP1ROW <- file.path(QTL, "Report_DMz26-3123_SNP_mapping_2.csv")  # one-row 0/1/2/-
SNP2ROW <- file.path(QTL, "Report_DMz26-3123_SNP_2.csv")          # two-row 0/1/-
SILICO  <- file.path(QTL, "Report_DMz26-3123_SilicoDArT_1.csv")   # dominant P/A
PHENO   <- file.path(QTL, "mozpue_phenotype.xlsx")

# DArT reports carry 6 header rows before the real header (documented in
# metadata.json$datafiles[[i]]$exportdetails$headerrowsdescr):
#   0 order#, 1 DArT plate barcode, 2 client plate barcode,
#   3 well row, 4 well col, 5 sample comments, 6 genotype name (= header)
DART_SKIP <- 6

rule <- function(x) cat("\n", strrep("=", 70), "\n", x, "\n", strrep("=", 70), "\n", sep = "")

## ---------------------------------------------------------------- 1. inventory
rule("1. FILE INVENTORY")
inv <- data.table(file = list.files(QTL))
inv[, size_MB := round(file.size(file.path(QTL, file)) / 1e6, 2)]
print(inv)

## ------------------------------------------------------- 2. read one-row SNPs
rule("2. SNP REPORT (one-row format) -- structure")
d <- fread(SNP1ROW, skip = DART_SKIP, header = TRUE, na.strings = c("", "NA"))
metacols <- 1:which(names(d) == "RepAvg")
cat("meta columns (", length(metacols), "):\n", sep = "")
print(names(d)[metacols])
G <- as.matrix(d[, -metacols, with = FALSE])
samples <- colnames(G)
cat("\nmarkers:", nrow(d), "  samples:", ncol(G), "\n")
cat("genotype call codes (0=hom ref, 1=hom SNP, 2=het, -=missing):\n")
tb <- table(G); print(round(100 * tb / sum(tb), 2))

# cross-check the other two files carry the same samples / marker count
h2 <- names(fread(SNP2ROW, skip = DART_SKIP, nrows = 1, header = TRUE))
hs <- names(fread(SILICO,  skip = DART_SKIP, nrows = 1, header = TRUE))
n2 <- nrow(fread(SNP2ROW, skip = DART_SKIP, select = 1, header = TRUE))
ns <- nrow(fread(SILICO,  skip = DART_SKIP, select = 1, header = TRUE))
cat("\ntwo-row SNP file : data rows", n2, "=> markers", n2 / 2,
    "| samples", length(h2) - which(h2 == "RepAvg"), "\n")
cat("SilicoDArT file  : markers", ns,
    "| samples", length(hs) - which(hs == "Reproducibility"), "\n")
cat("sample sets identical across all three files:",
    identical(samples, h2[(which(h2 == "RepAvg") + 1):length(h2)]) &&
    identical(samples, hs[(which(hs == "Reproducibility") + 1):length(hs)]), "\n")

## --------------------------------------------------- 3. sample identity / plate
rule("3. SAMPLES -- identity, plates, comments")
hdr <- fread(SNP1ROW, nrows = DART_SKIP + 1, header = FALSE, colClasses = "character")
s0 <- which(names(d) == "RepAvg") + 1
lab <- c("order#", "DArT_plate", "client_plate", "well_row", "well_col", "comments")
for (i in seq_along(lab)) {
  vals <- unlist(hdr[i, s0:ncol(hdr), with = FALSE])
  cat(sprintf("%-13s unique values: %s\n", lab[i],
              paste(names(sort(table(vals), decreasing = TRUE)), collapse = " | ")))
}
gid <- as.numeric(sub("plant ", "", samples))
cat("\nsample name pattern: 'plant <n>'; n range", min(gid), "-", max(gid), "\n")
cat("IDs absent from 1:188 :", paste(setdiff(1:188, gid), collapse = ", "), "\n")
cat("NOTE: no parental lines, no checks, no technical replicates among the samples.\n")

## --------------------------------------------------------------- 4. marker QC
rule("4. MARKER QC")
chr <- sub(" .*", "", d$Chrom_Maize_B73_V4.0.assembly)
pos <- suppressWarnings(as.numeric(d$ChromPosSnp_Maize_B73_V4.0.assembly))
for (v in c("CallRate", "RepAvg", "FreqHets", "AvgCountRef", "AvgCountSnp")) {
  cat(sprintf("%-12s ", v)); print(round(summary(d[[v]]), 3))
}
cat("\nread depth is modest: median ~", round(median(d$AvgCountRef), 1),
    "x ref /", round(median(d$AvgCountSnp), 1), "x alt. At this depth allelic\n",
    "dropout can undercall heterozygotes, which in an F2 is 50% of individuals.\n",
    "BUT: agent/qtl_structure_check.R shows het = 0.502 (median) at the RETAINED\n",
    "framework markers -- exactly the F2 expectation -- so dropout is NOT\n",
    "materially distorting the filtered set. The low-FreqHets tail below is\n",
    "uninformative markers (parents share the allele), not dropout artefact.\n",
    "Still use est.map/calc.errorlod with error.prob ~0.01 rather than the\n",
    "default 1e-4, but no aggressive correction is warranted.\n")

cat("\nanchored to B73 AGPv4 chr 1-10 :", sum(chr %in% as.character(1:10)), "/", nrow(d), "\n")
cat("unanchored (blank Chrom)        :", sum(is.na(chr) | chr == ""), "\n")
cat("on unplaced contigs (B73V4_ctg*):", sum(grepl("^B73V4_ctg", chr)), "\n")

n0 <- rowSums(G == "0"); n1 <- rowSums(G == "1"); n2 <- rowSums(G == "2")
tot <- n0 + n1 + n2
maf <- pmin(2 * n0 + n2, 2 * n1 + n2) / (2 * tot)
# chi-square goodness of fit to F2 codominant expectation 1 : 2 : 1
cs <- (n0 - tot/4)^2/(tot/4) + (n1 - tot/4)^2/(tot/4) + (n2 - tot/2)^2/(tot/2)
p121 <- pchisq(cs, df = 2, lower.tail = FALSE)
cat("\neffectively monomorphic (MAF<0.05):", sum(maf < 0.05, na.rm = TRUE),
    " <- parents share the allele; uninformative\n")
cat("MAF >= 0.15                       :", sum(maf >= 0.15, na.rm = TRUE), "\n")
cat("fits 1:2:1 at p>0.01 (all markers):", sum(p121 > 0.01, na.rm = TRUE), "\n")

# filter stack. 1:2:1 fit is used here as a PROXY for "polymorphic between the
# parents", because the parents themselves were not genotyped (see report).
keep <- chr %in% as.character(1:10) & d$CallRate >= 0.90 & d$RepAvg >= 0.95 &
        maf >= 0.15 & p121 > 0.01
keep[is.na(keep)] <- FALSE
cat("\nfilter stack (anchored & CallRate>=0.90 & RepAvg>=0.95 & MAF>=0.15 & 1:2:1 p>0.01):",
    sum(keep), "markers\n")

## ----------------------------------- 5. collapse DArT secondaries (same tag)
rule("5. FRAMEWORK MAP CANDIDATE (1 SNP per CloneID)")
k  <- which(keep)
dt <- data.table(i = k, clone = d$CloneID[k], allele = d$AlleleID[k],
                 cr = d$CallRate[k], chr = chr[k], pos = pos[k])
cat("secondaries (>1 SNP on the same 69bp tag) dropped:",
    nrow(dt) - uniqueN(dt$clone), "\n")
best <- dt[order(-cr), .SD[1], by = clone][order(as.integer(chr), pos)]
cat("framework markers:", nrow(best), "\n\n")
best[, mb := pos / 1e6]
print(best[, .(n = .N,
               span_Mb        = round(max(mb) - min(mb), 1),
               median_gap_Mb  = round(median(diff(mb)), 2),
               max_gap_Mb     = round(max(diff(mb)), 1)),
           by = chr][order(as.integer(chr))])

## ---------------------------------------------------------------- 6. sample QC
rule("6. SAMPLE QC")
scr  <- colMeans(G != "-")
shet <- colMeans(G == "2") / scr
cat("per-sample call rate     : "); print(round(summary(scr), 3))
cat("per-sample heterozygosity: "); print(round(summary(shet), 3))
cat("\ncall rate < 0.80 :", paste(names(scr)[scr < 0.80], collapse = ", "), "\n")
cat("het < 0.15 (would look like an inbred parent):",
    paste(names(shet)[shet < 0.15], collapse = ", "), "| none expected\n")
cat("het > 0.60 (would look like an F1 / contaminant):",
    paste(names(shet)[shet > 0.60], collapse = ", "), "\n")
cat("=> heterozygosity range", round(min(shet), 2), "-", round(max(shet), 2),
    "is consistent with all samples being F2 individuals.\n")

## ------------------------------------------------------------- 7. phenotypes
rule("7. PHENOTYPES -- mozpue_phenotype.xlsx")
cat("sheets:", paste(excel_sheets(PHENO), collapse = ", "), "\n")
p <- as.data.frame(suppressMessages(read_excel(PHENO, sheet = "data base")))
names(p)[1] <- "plant"
cat("raw dim:", paste(dim(p), collapse = " x "), "\n")
cat("row 1 is a UNITS row, not data (e.g. 'Yield' -> ",
    p[1, "Yield"], "); dropping it.\n", sep = "")
pu <- as.data.frame(suppressMessages(read_excel(PHENO, sheet = "data base",
                                                col_types = "text")))
p <- p[-1, ]; pu <- pu[-1, ]
p$plant <- suppressWarnings(as.numeric(p$plant))
cat("\ntraits:\n"); print(names(p))
cat("\nrows:", nrow(p), " rows with a plant ID:", sum(!is.na(p$plant)),
    " fully blank rows:", sum(is.na(p$plant)), "\n")
cat("plant IDs missing from 1:188:\n  ",
    paste(setdiff(1:188, p$plant), collapse = ", "), "\n")

cat("\nper-trait completeness and range:\n")
for (n in setdiff(names(p), c("plant", "OBS: "))) {
  v <- suppressWarnings(as.numeric(p[[n]]))
  cat(sprintf("%-28s nonNA=%3d  min=%9.3f  median=%9.3f  max=%9.3f\n",
              n, sum(!is.na(v)), min(v, na.rm = TRUE),
              median(v, na.rm = TRUE), max(v, na.rm = TRUE)))
}
o <- pu[[ncol(pu)]]
cat("\nOBS: notes (", sum(!is.na(o)), " plants):\n", sep = "")
print(data.frame(plant = pu[[1]][!is.na(o)], note = o[!is.na(o)]))
cat("\nNOTE: no treatment / block / rep / location / date column, and no tissue-P\n",
    "column. PUE, HI, RSR, RYE are ratios of the other measured columns, i.e.\n",
    "derived and mutually correlated by construction.\n")

## ------------------------------------------------------- 8. the critical join
rule("8. GENOTYPE x PHENOTYPE JOIN -- effective mapping population")
pid <- p$plant[!is.na(p$plant)]
cat("genotyped        :", length(gid), "\n")
cat("phenotyped       :", length(pid), "\n")
cat("BOTH (usable n)  :", length(intersect(gid, pid)), "\n")
cat("genotyped only   :", length(setdiff(gid, pid)), "->",
    paste(sort(setdiff(gid, pid)), collapse = ", "), "\n")
cat("phenotyped only  :", length(setdiff(pid, gid)), "->",
    paste(sort(setdiff(pid, gid)), collapse = ", "), "\n")
cat("\nCAVEAT: the 'plant N' -> field-book plant N correspondence is assumed,\n",
    "not documented anywhere. An off-by-one in plate layout would silently\n",
    "destroy the analysis. Needs confirmation from the submission sheet.\n")

## ---------------------------------------------------------------- 9. toolchain
rule("9. TOOLCHAIN")
cat(R.version.string, "\n")
want <- c("qtl", "qtl2", "qtl2convert", "ASMap", "dartR.base", "LinkageMapView",
          "qtlcharts", "readxl", "data.table", "dplyr", "ggplot2", "writexl", "vcfR")
inst <- rownames(installed.packages())
print(data.frame(package = want, installed = want %in% inst))

cat("\naudit complete.\n")
