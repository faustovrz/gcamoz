#!/usr/bin/env Rscript
# qtl_make_abh.R -- convert the DArTseq F2 calls to A/H/B for R/qtl.
#
# PROVENANCE: reimplements TASSEL 5 GenosToABHPlugin. Algorithm read from
#   https://raw.githubusercontent.com/zacharymiller90/tassel-ML/refs/heads/master/
#     src/net/maizegenetics/analysis/data/GenosToABHPlugin.java   (md5 2f72573e4a…)
# Retained ONLY if agent/qtl_abh_parity_test.sh shows parity against the local
# /Applications/TASSEL 5/run_pipeline.pl binary (which differs from that fork --
# it accepts -outputFormat, the fork has no such parameter).
#
# ARCHITECTURE: abh_call() below operates on diploid genotype STRINGS ("A/A",
# "A/G", "N/N") -- the same abstraction TASSEL works on -- so the parity harness
# exercises exactly the function that runs on real data. The DArT path is a thin
# adapter that builds those strings from 0/1/2/- codes plus the CML530 allele.
#
# FIDELITY TO TASSEL (GenosToABHPlugin.java line refs):
#   L157-166  site rejected if: pA==pB | pA unknown | pB unknown | pA het | pB het
#   L173-174  H matches EITHER phase of the parental heterozygote
#   L249-258  exact diploid equality, tested in order A, B, H, else NA
#   L191-212  parental consensus is UNANIMOUS-OR-UNKNOWN (any conflict -> unknown).
#             Not a modal vote. Irrelevant here (one CML530 observation, no
#             replicates) but implemented for the parity test.
#
# DELIBERATE DIVERGENCES (justified, and reported so they are auditable):
#   1. TASSEL collapses missing data AND off-parent third alleles both into "NA"
#      (the fall-through at L257). We have that information, so they are tracked
#      separately: missing -> NA, off-parent -> counted and reported.
#   2. TASSEL DROPS rejected sites from the output entirely (L222-224). We keep
#      unpolarizable markers with ARBITRARY orientation, flagged in
#      polarity_source. Recombination fractions are orientation-invariant, so
#      such markers still carry full mapping information; only the sign of the
#      allele effect is unavailable. Dropping them would discard real data.
#
# Runs BEFORE the CML530 genotypes exist: if CML530_marker_alleles.tsv is absent
# every marker gets polarity_source="arbitrary", so the R/qtl format and the
# read.cross() round-trip can be validated now and re-run for polarity later.
#
# Usage: Rscript agent/qtl_make_abh.R > agent/qtl_make_abh.log 2>&1

suppressMessages({library(data.table); library(readxl); library(qtl)})

SNP1ROW   <- "data/qtl/Report_DMz26-3123_SNP_mapping_2.csv"
DART_SKIP <- 6
DERIVED   <- "data/qtl/derived"
ROSTER    <- file.path(DERIVED, "markers_v5.tsv")
# Prefer the DIRECT assembly read (qtl_cml530_alleles_direct.R) over the pileup
# route: every pileup site was ~DP=1, where bcftools call -m defaults toward REF
# and produced random polarity (53% phase-consistent, agent/qtl_diagnose_phase.log).
CMLCALLS  <- if (file.exists(file.path(DERIVED, "CML530_marker_alleles_direct.tsv")))
               file.path(DERIVED, "CML530_marker_alleles_direct.tsv") else
               file.path(DERIVED, "CML530_marker_alleles.tsv")
PHENO     <- "data/qtl/mozpue_phenotype.xlsx"
DICT      <- "data/phenotype_dictionary.csv"   # single authority for trait names
GENFILE   <- file.path(DERIVED, "rqtl_gen_abh.csv")
PHEFILE   <- file.path(DERIVED, "rqtl_phe.csv")
MKINFO    <- file.path(DERIVED, "abh_marker_info.tsv")

rule <- function(x) cat("\n", strrep("=", 70), "\n", x, "\n", strrep("=", 70), "\n", sep = "")

## ===================================================== TASSEL-equivalent core
UNK <- "N/N"

#' Normalise a diploid genotype string to an unordered allele pair.
.norm <- function(g) {
  g[is.na(g) | g == "" | g == "./." | g == "NA"] <- UNK
  a <- toupper(sub("[/|].*$", "", g)); b <- toupper(sub("^.*[/|]", "", g))
  a[a %in% c("N", ".", "")] <- "N"; b[b %in% c("N", ".", "")] <- "N"
  ifelse(a == "N" | b == "N", UNK, paste0(pmin(a, b), "/", pmax(a, b)))
}
.is_het <- function(g) { g <- .norm(g); g != UNK & sub("/.*", "", g) != sub(".*/", "", g) }
.is_unk <- function(g) .norm(g) == UNK

#' TASSEL GenosToABHPlugin site acceptance test (L157-166).
#' @return TRUE if the site is usable for ABH coding.
abh_site_ok <- function(pA, pB) {
  pA <- .norm(pA); pB <- .norm(pB)
  !(pA == pB | .is_unk(pA) | .is_unk(pB) | .is_het(pA) | .is_het(pB))
}

#' TASSEL parental consensus (L191-212): unanimous-or-unknown, NOT a modal vote.
#' Missing calls are skipped; two conflicting non-missing calls -> unknown.
abh_consensus <- function(gts) {
  g <- .norm(gts); g <- g[!is.na(g)]
  if (!length(g)) return(UNK)
  fin <- g[1]; unk <- fin == UNK
  for (i in seq_along(g)[-1]) {
    if (g[i] != fin) {
      if (unk) { fin <- g[i]; unk <- FALSE }
      else if (g[i] != UNK) return(UNK)   # conflict aborts immediately
    }
  }
  fin
}

#' TASSEL genotype -> ABH call (L249-258). Exact diploid equality, ordered.
#' Returns "A", "B", "H", or NA. `off` receives TRUE where the genotype is a
#' real call that matches none of A/B/H (TASSEL merges these into NA; we do not).
abh_call <- function(geno, pA, pB) {
  g <- .norm(geno); a <- .norm(pA); b <- .norm(pB)
  het <- .norm(paste0(sub("/.*", "", a), "/", sub("/.*", "", b)))  # both phases
  out <- rep(NA_character_, length(g))
  out[g == a]   <- "A"
  out[g == b]   <- "B"
  out[g == het] <- "H"
  attr(out, "off_parent") <- is.na(out) & g != UNK
  out
}

## ------------------------------------------------------------------ 1. inputs
rule("1. INPUTS")
d <- fread(SNP1ROW, skip = DART_SKIP, header = TRUE, na.strings = c("", "NA"))
mc <- 1:which(names(d) == "RepAvg")
G  <- as.matrix(d[, -mc, with = FALSE]); rownames(G) <- d$AlleleID
roster <- fread(ROSTER, colClasses = list(character = "chr_v5"))
cat("DArT markers:", nrow(d), " samples:", ncol(G), "\n")
cat("roster (lifted to v5):", nrow(roster), " framework:", sum(roster$framework), "\n")

have_cml <- file.exists(CMLCALLS)
if (have_cml) {
  cml <- fread(CMLCALLS, colClasses = list(character = "chr_v5"))
  cat("CML530 calls:", nrow(cml), "\n"); print(cml[, .N, by = cml_class][order(-N)])
} else {
  cat("\n*** CML530_marker_alleles.tsv NOT FOUND ***\n")
  cat("Running in ARBITRARY-POLARITY mode: every marker keeps a consistent but\n")
  cat("unattributed orientation. Map/QTL positions are unaffected; allele-effect\n")
  cat("SIGNS are not interpretable. Re-run after the CML530 pipeline finishes.\n")
  cml <- data.table(marker = character(), cml_class = character())
}

## --------------------------------------- 2. DArT codes -> genotype strings
rule("2. ADAPTER: DArT 0/1/2/- -> diploid genotype strings")
# DArT `SNP` column is e.g. "28:A>G": offset, then REF>ALT. REF is the B73 allele
# (DArT called it against AGPv4), NOT a statement about either parent.
mk <- roster[marker %in% rownames(G)]
mk[, ref_allele := sub(".*:([ACGT])>([ACGT]).*", "\\1", SNP)]
mk[, alt_allele := sub(".*:([ACGT])>([ACGT]).*", "\\2", SNP)]
bad <- mk[!(ref_allele %in% c("A","C","G","T")) | !(alt_allele %in% c("A","C","G","T"))]
cat("markers with unparseable SNP field:", nrow(bad), "\n")
mk <- mk[ref_allele %in% c("A","C","G","T") & alt_allele %in% c("A","C","G","T")]

mk <- merge(mk, cml[, .(marker, cml_class)], by = "marker", all.x = TRUE)
mk[is.na(cml_class), cml_class := "no-call"]
# Polarity: parent A := the CML530 allele. Where CML530 is unresolved (no-hit,
# het-call = paralog signal, or not yet computed) fall back to ref as parent A
# and record that the orientation is arbitrary.
mk[, polarity_source := fifelse(cml_class %in% c("ref", "alt"), "CML530", "arbitrary")]
mk[, pA := fifelse(cml_class == "alt", alt_allele, ref_allele)]
mk[, pB := fifelse(cml_class == "alt", ref_allele, alt_allele)]
mk[, `:=`(pA_gt = paste0(pA, "/", pA), pB_gt = paste0(pB, "/", pB))]
print(mk[, .N, by = .(cml_class, polarity_source)][order(-N)])

## ------------------------------------------------------------ 3. ABH encoding
rule("3. ABH ENCODING (TASSEL-equivalent core)")
mk <- mk[framework == TRUE]
setorder(mk, chr_v5, pos_v5)
cat("framework markers encoded:", nrow(mk), "\n")
Gm <- G[mk$marker, , drop = FALSE]
# DArT one-row: 0 = hom REF, 1 = hom ALT, 2 = het, - = missing
gt <- matrix(UNK, nrow(Gm), ncol(Gm), dimnames = dimnames(Gm))
gt[Gm == "0"] <- paste0(mk$ref_allele, "/", mk$ref_allele)[row(Gm)[Gm == "0"]]
gt[Gm == "1"] <- paste0(mk$alt_allele, "/", mk$alt_allele)[row(Gm)[Gm == "1"]]
gt[Gm == "2"] <- paste0(mk$ref_allele, "/", mk$alt_allele)[row(Gm)[Gm == "2"]]

site_ok <- abh_site_ok(mk$pA_gt, mk$pB_gt)
cat("sites passing the TASSEL acceptance test:", sum(site_ok), "of", nrow(mk), "\n")
cat("  (all should pass: parents are constructed homozygous and distinct)\n")

ABH <- matrix(NA_character_, nrow(Gm), ncol(Gm), dimnames = dimnames(Gm))
off <- matrix(FALSE, nrow(Gm), ncol(Gm))
for (i in seq_len(nrow(Gm))) {
  cl <- abh_call(gt[i, ], mk$pA_gt[i], mk$pB_gt[i])
  ABH[i, ] <- cl; off[i, ] <- attr(cl, "off_parent")
}
tb <- table(factor(ABH, levels = c("A", "H", "B")), useNA = "no")
cat("\nABH composition:\n"); print(tb)
cat(sprintf("  A %.1f%%  H %.1f%%  B %.1f%%   (F2 expectation 25/50/25)\n",
            100*tb["A"]/sum(tb), 100*tb["H"]/sum(tb), 100*tb["B"]/sum(tb)))
cat("missing (NA):", sum(is.na(ABH)), sprintf("(%.1f%%)\n", 100*mean(is.na(ABH))))
cat("off-parent calls (TASSEL would silently merge these into NA):", sum(off), "\n")

## -------------------------------------------------------- 4. R/qtl csvsr out
rule("4. WRITE R/qtl csvsr FILES")
p <- as.data.frame(suppressMessages(read_excel(PHENO, sheet = "data base")))
names(p)[1] <- "plant"
names(p) <- trimws(names(p))   # the "OBS:" header carries stray trailing space
p <- p[-1, ]                   # row 1 is a units row, not data
p$plant <- suppressWarnings(as.numeric(p$plant))
p <- p[!is.na(p$plant), ]

# Trait names MUST be the publication abbreviations, per the convention set in
# male_gca_trait_sweep.qmd. data/phenotype_dictionary.csv is the single authority
# (built by agent/phenotype_dictionary_unify.R); its raw_header column maps each
# abbreviation to the literal Excel column, so the mapping is NOT duplicated here.
dict <- fread(DICT)[project == "qtl_f2"]
tr <- dict[role == "trait"]
missing_hdr <- setdiff(tr$raw_header, names(p))
if (length(missing_hdr)) stop("dictionary raw_header not found in phenotype file: ",
                              paste(missing_hdr, collapse = " | "))
setnames(p, tr$raw_header, tr$variable)
cat("traits renamed to abbreviations:", paste(tr$variable, collapse = " "), "\n")
unmapped <- setdiff(names(p), c("plant", tr$variable, dict[role != "trait", variable],
                                dict[role != "trait", raw_header]))
if (length(unmapped))
  cat("*** columns absent from the dictionary (ADD THEM):",
      paste(unmapped, collapse = " | "), "***\n")

gid <- as.numeric(sub("plant ", "", colnames(ABH)))
keep_ind <- which(gid %in% p$plant)             # 166 with both geno and pheno
ids <- paste0("plant", gid[keep_ind])
cat("individuals with genotype AND phenotype:", length(keep_ind), "\n")

# csvsr genotype file: one row per marker -> id, chr, then one col per individual.
# EXACTLY TWO leading columns. Verified empirically by agent/qtl_rqtl_format_probe.R:
# a third `pos` column is parsed as an EXTRA INDIVIDUAL named "pos" (giving 167
# instead of 166) and its values become unknown genotype codes -- a silent
# corruption, not an error. csvsr carries no positions; R/qtl assigns dummy 5 cM
# spacing, which est.map replaces. Physical positions live in MKINFO instead.
gen <- data.table(id = mk$marker, chr = mk$chr_v5)
gen <- cbind(gen, as.data.table(ABH[, keep_ind, drop = FALSE]))
setnames(gen, c("id", "chr", ids))
fwrite(gen, GENFILE, na = "-", quote = FALSE)

# csvsr phenotype file is ALSO rotated: one row per phenotype. Trait ORDER comes
# from the dictionary so it is stable across runs rather than following Excel.
traits <- tr$variable
pm <- p[match(gid[keep_ind], p$plant), traits, drop = FALSE]
phe <- data.table(id = traits)
phe <- cbind(phe, as.data.table(t(sapply(traits, function(t)
  suppressWarnings(as.numeric(pm[[t]]))))))
setnames(phe, c("id", ids))
fwrite(phe, PHEFILE, na = "-", quote = FALSE)
cat("wrote", GENFILE, "  dim:", paste(dim(gen), collapse = " x "), "\n")
cat("wrote", PHEFILE, "  dim:", paste(dim(phe), collapse = " x "), "\n")

fwrite(mk[, .(marker, CloneID, chr_v5, pos_v5, chr_v4, pos_v4, SNP,
              ref_allele, alt_allele, cml_class, pA, pB, polarity_source,
              CallRate, RepAvg, maf, p121, framework)], MKINFO, sep = "\t")
cat("wrote", MKINFO, "\n")

## ------------------------------------------------- 5. validate via read.cross
rule("5. VALIDATE: read.cross() round-trip")
# The real format test: if R/qtl can load it, the layout is right. Assert rather
# than claim.
cr <- read.cross(format = "csvsr", dir = "", genfile = GENFILE, phefile = PHEFILE,
                 genotypes = c("A", "H", "B"), na.strings = c("-", "NA"),
                 crosstype = "f2", estimate.map = FALSE)
print(summary(cr))
cat("\nnind:", nind(cr), " ntotmar:", totmar(cr), " nchr:", nchr(cr), "\n")
cat("nphe:", nphe(cr), "\n")
stopifnot(nind(cr) == length(keep_ind), totmar(cr) == nrow(mk), nchr(cr) == 10)
cat("\nread.cross() round-trip OK\n")

rule("NOTES")
cat("1. `pos` is v5 Mb, a PLACEHOLDER. Genetic distances come from est.map on\n")
cat("   the physical order, per the agreed approach.\n")
cat("2. polarity_source in", basename(MKINFO), "records per marker whether A is\n")
cat("   the CML530 allele or an arbitrary orientation. The assumption travels\n")
cat("   with the data rather than living in a comment.\n")
cat("3. NOT YET VALIDATED against TASSEL. Run agent/qtl_abh_parity_test.sh and\n")
cat("   discard this script if it does not match /Applications/TASSEL 5/.\n")
