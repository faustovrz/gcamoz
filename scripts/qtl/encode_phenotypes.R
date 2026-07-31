#!/usr/bin/env Rscript
# encode_phenotypes.R -- write the R/qtl csvsr PHENOTYPE file.
#
# Split out of make_abh_framework_only.R (the superseded 1,384-marker encoder, which is
# on the delete list). The phenotype writer had no business living inside a genotype
# encoder, and it meant rqtl_phe.csv was produced by a script slated for deletion.
#
# TWO CHANGES from the version it replaces:
#   1. Individuals are taken from the GENOTYPE file, so the two files always carry the
#      same individuals in the same order. The old version emitted only the 166 with
#      both genotype and phenotype, while the genotype file has 186 -- read.cross then
#      warned about 20 mismatches on every load.
#   2. The 20 genotyped-but-unphenotyped plants get NA columns rather than being
#      dropped. est.map does not use phenotypes, so the map should use all 186
#      genotyped plants; scanone will drop the NAs per trait by itself.
#
# Trait NAMES come from data/phenotype_dictionary.csv (the single authority, built by
# phenotype_dictionary_unify.R) via its raw_header column, so the Excel-header ->
# abbreviation mapping is not duplicated here.
#
# csvsr phenotype layout is ROTATED: one row per phenotype, one column per individual.
#
# Reads   data/qtl/mozpue_phenotype.xlsx
#         data/phenotype_dictionary.csv
#         data/qtl/derived/rqtl_gen_abh_all.csv   (for the individual roster only)
# Writes  data/qtl/derived/rqtl_phe.csv
# Usage: Rscript scripts/qtl/encode_phenotypes.R > agent/encode_phenotypes.log 2>&1

suppressMessages({library(data.table); library(readxl); library(qtl)})
rule <- function(x) cat("\n", strrep("=", 70), "\n", x, "\n", strrep("=", 70), "\n", sep = "")

D       <- "data/qtl/derived"
PHENO   <- "data/qtl/mozpue_phenotype.xlsx"
DICT    <- "data/phenotype_dictionary.csv"
GENFILE <- file.path(D, "rqtl_gen_abh_all.csv")
OUT     <- file.path(D, "rqtl_phe.csv")

rule("1. INDIVIDUAL ROSTER (from the genotype file, so the two always match)")
stopifnot(file.exists(GENFILE))
ids <- names(fread(GENFILE, sep = ",", nrows = 0))[-(1:2)]   # drop id, chr
cat("individuals in the genotype file:", length(ids), "\n")
cat("first five:", paste(head(ids, 5), collapse = " "), "\n")
gid <- suppressWarnings(as.numeric(sub("^plant", "", ids)))
stopifnot(!any(is.na(gid)))

rule("2. PHENOTYPES")
p <- as.data.frame(suppressMessages(read_excel(PHENO, sheet = "data base")))
names(p)[1] <- "plant"
names(p) <- trimws(names(p))     # the "OBS:" header carries stray trailing space
p <- p[-1, ]                     # row 1 is a UNITS row, not data
p$plant <- suppressWarnings(as.numeric(p$plant))
p <- p[!is.na(p$plant), ]
cat("rows with a plant ID:", nrow(p), "\n")

dict <- fread(DICT)[project == "qtl_f2"]
tr   <- dict[role == "trait"]
missing_hdr <- setdiff(tr$raw_header, names(p))
if (length(missing_hdr))
  stop("dictionary raw_header not found in the phenotype file: ",
       paste(missing_hdr, collapse = " | "))
setnames(p, tr$raw_header, tr$variable)
cat("traits (from the dictionary):", paste(tr$variable, collapse = " "), "\n")
unmapped <- setdiff(names(p), c("plant", tr$variable,
                                dict[role != "trait", variable],
                                dict[role != "trait", raw_header]))
if (length(unmapped))
  cat("*** columns absent from the dictionary (ADD THEM):",
      paste(unmapped, collapse = " | "), "***\n")

rule("3. ALIGN TO THE GENOTYPE ROSTER")
idx <- match(gid, p$plant)                # NA where a genotyped plant has no phenotype
cat("genotyped individuals            :", length(gid), "\n")
cat("with phenotypes                  :", sum(!is.na(idx)), "\n")
cat("genotyped but NOT phenotyped     :", sum(is.na(idx)), "\n")
if (any(is.na(idx)))
  cat("  ->", paste(sort(gid[is.na(idx)]), collapse = ", "), "\n")
cat("phenotyped but NOT genotyped     :", sum(!p$plant %in% gid), "\n")
if (any(!p$plant %in% gid))
  cat("  ->", paste(sort(p$plant[!p$plant %in% gid]), collapse = ", "), "\n")

rule("4. WRITE csvsr (rotated: one ROW per phenotype)")
mat <- t(sapply(tr$variable, function(v) suppressWarnings(as.numeric(p[[v]]))[idx]))
phe <- cbind(data.table(id = tr$variable), as.data.table(mat))
setnames(phe, c("id", ids))
stopifnot(nrow(phe) == nrow(tr), ncol(phe) == length(ids) + 1L)
fwrite(phe, OUT, na = "-", quote = FALSE)
cat("wrote", OUT, " dim:", paste(dim(phe), collapse = " x "), "\n")
cat("\nper-trait completeness:\n")
print(data.table(trait = tr$variable,
                 n_nonNA = apply(mat, 1, function(x) sum(is.finite(x))),
                 median  = round(apply(mat, 1, median, na.rm = TRUE), 3)))

rule("5. VERIFY read.cross ROUND-TRIP")
cr <- read.cross(format = "csvsr", dir = "", genfile = GENFILE, phefile = OUT,
                 genotypes = c("A","H","B"), na.strings = c("-","NA"),
                 crosstype = "f2", estimate.map = FALSE)
cat("nind:", nind(cr), " ntotmar:", totmar(cr), " nchr:", nchr(cr),
    " nphe:", nphe(cr), "\n")
stopifnot(nind(cr) == length(ids))
cat("\nindividual sets match the genotype file -- no read.cross mismatch warning.\n")
