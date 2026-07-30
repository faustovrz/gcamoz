#!/usr/bin/env bash
# qtl_cml530_genotypes.sh -- CML530 genotypes at the 3,828 DArTseq marker
# positions, in B73 v5 coordinates.
#
# DESIGN (per FVRZ, 2026-07-30):
#   The join key is the B73 v5 COORDINATE, not a read name. The markers were
#   already lifted to v5 (agent/qtl_liftover_v4_to_v5.R -> markers_v5.tsv), so
#   the final call is a plain pileup at known sites. Read names are irrelevant.
#
#   Two mapping steps, because the CML530 sequence has to be obtained before it
#   can be placed in v5:
#     [1] tags -> CML530 assembly     : locate each tag's homologous locus
#     [2] CML530 span -> B73 v5       : put CML530's own sequence in v5 coords
#     [3] pileup at markers_v5.bed    : read the genotype positionally
#
#   Extraction in step [2] is the ALIGNED SPAN ONLY -- no flanking sequence. The
#   tag is the assay; extra flank adds nothing to the call. Paralogous/ambiguous
#   placement is detected from MAPQ / NM / secondary hits in the BAM, not by
#   re-mapping a longer window.
#
#   Genotype is taken from the pileup rather than by walking the SNP offset by
#   hand: if CML530 carries an indel near the site, the CIGAR keeps the base
#   assignment honest where offset arithmetic would silently shift it.
#
# CML530 IS A PARENT OF THIS POPULATION -- not a proxy. Confirmed from
# README.md:44-46 and multilocation_alpha_lattice.qmd:62, which list the NC II
# female set as CML364, CML366, CML434, CML435, CML439, CML530, CML532, and note
# explicitly that "the earlier draft had incorrectly listed CML440". Stelio's
# July 2025 email was that earlier draft; the repo supersedes it. So the allele
# read here is the ACTUAL parental allele and ABH polarity is directly
# determined, not proxied.
#
# GENERALISES: one BAM per accession. The same tag FASTA works unchanged for the
# ENxx tester (unsequenced) or any other CML female, should those be sequenced.
#
# Usage:
#   bash agent/qtl_cml530_genotypes.sh --preflight     # check inputs only, no compute
#   bash agent/qtl_cml530_genotypes.sh                 # full run
# Log:
#   bash agent/qtl_cml530_genotypes.sh > agent/qtl_cml530_genotypes.log 2>&1

set -euo pipefail

# ---- configuration -----------------------------------------------------------
REPO="/Users/fvrodriguez/repos/gcamoz"
ZEAL="/Users/fvrodriguez/repos/zealhmm"

ACCESSION="CML530"
CML_GZ="${ZEAL}/data/ref/genomes/Zm-CML530-REFERENCE-HiLo-1.0/Zm-CML530-REFERENCE-HiLo-1.0.fa.gz"
# CML530 ships as BGZF (block-compressed, gzip-compatible), so samtools/minibwa
# both read it in place. NO decompression step -- saves 2.3 GB of disk.
# B73 NAM v5 already on disk (found 2026-07-30); reused read-only, not re-downloaded.
B73V5_FA="/Users/fvrodriguez/repos/zealtiger/data/brbseq/ref/Zm-B73-REFERENCE-NAM-5.0.fa"

# minibwa index prefixes live here so indexing never writes into the source ref
# directories (B73 v5 belongs to zealtiger; CML530 to zealhmm).
IDXDIR="${ZEAL}/data/ref/genomes/minibwa_index"
CML_IDX="${IDXDIR}/${ACCESSION}"
B73_IDX="${IDXDIR}/B73v5"

DERIVED="${REPO}/data/qtl/derived"
TAGS="${DERIVED}/tags_for_cml530.fa"        # from qtl_liftover_v4_to_v5.R
ROSTER="${DERIVED}/markers_v5.tsv"          # from qtl_liftover_v4_to_v5.R
BED="${DERIVED}/markers_v5.bed"             # written by step 0 below
OUT="${DERIVED}/${ACCESSION}"
THREADS="${THREADS:-8}"

PREFLIGHT_ONLY=0
[[ "${1:-}" == "--preflight" ]] && PREFLIGHT_ONLY=1

rule() { printf '\n%s\n%s\n%s\n' "======================================================================" "$1" "======================================================================"; }

# ---- 0. preflight ------------------------------------------------------------
rule "0. PREFLIGHT"
fail=0
for t in minibwa samtools bcftools Rscript; do
  if command -v "$t" >/dev/null 2>&1; then
    printf 'tool  OK      %-10s %s\n' "$t" "$(command -v "$t")"
  else
    printf 'tool  MISSING %-10s\n' "$t"; fail=1
  fi
done
for f in "$CML_GZ" "$TAGS" "$ROSTER"; do
  if [[ -s "$f" ]]; then
    printf 'input OK      %8s  %s\n' "$(du -h "$f" | cut -f1)" "$f"
  else
    printf 'input MISSING           %s\n' "$f"; fail=1
  fi
done
if [[ -s "$B73V5_FA" ]]; then
  printf 'input OK      %8s  %s\n' "$(du -h "$B73V5_FA" | cut -f1)" "$B73V5_FA"
else
  printf 'input MISSING           %s\n' "$B73V5_FA"
  printf '      ^ B73 NAM v5 FASTA needed for step 2. Not yet on disk.\n'
  fail=1
fi

# memory note: minibwa index uses 18N RAM (N = genome size). Maize ~2.2 Gb ->
# ~40 GB, and this box has 24 GB, so -l (low-memory GPL BWT) is MANDATORY here.
RAM_GB=$(( $(sysctl -n hw.memsize) / 1024 / 1024 / 1024 ))
printf '\nRAM: %s GB -> using "minibwa index -l" (18N would need ~40 GB)\n' "$RAM_GB"
printf 'threads: %s\n' "$THREADS"

if [[ "$fail" -ne 0 ]]; then
  printf '\nPREFLIGHT FAILED - resolve the MISSING items above before running.\n'
  exit 1
fi
printf '\npreflight OK\n'
[[ "$PREFLIGHT_ONLY" -eq 1 ]] && { printf 'stopping (--preflight)\n'; exit 0; }

mkdir -p "$OUT"

# ---- 1. marker BED in v5 coordinates ----------------------------------------
rule "1. MARKER BED (B73 v5 sites)"
# CRITICAL: markers_v5.tsv carries chromosomes as "1".."10", because that is what
# the v4->v5 chain file uses on its target side. BOTH FASTAs on disk name them
# "chr1".."chr10". Without the prefix, `bcftools mpileup -T` matches nothing and
# returns an EMPTY result rather than an error -- a silent failure. Verified
# 2026-07-30 from the .fai of each reference.
Rscript -e '
  suppressMessages(library(data.table))
  a <- commandArgs(TRUE); r <- fread(a[1], colClasses = list(character = "chr_v5"))
  bed <- r[, .(chrom = paste0("chr", chr_v5), start = pos_v5 - 1L, end = pos_v5,
               name = marker, score = round(100 * CallRate), strand = ".")]
  setorder(bed, chrom, start)
  fwrite(bed, a[2], sep = "\t", col.names = FALSE)
  cat("sites:", nrow(bed), " framework:", sum(r$framework), "\n")
  cat("chrom values written:", paste(sort(unique(bed$chrom)), collapse = " "), "\n")
' "$ROSTER" "$BED"
printf 'wrote %s\n' "$BED"

# ---- 2. CML530 reference (BGZF, indexed in place) ---------------------------
rule "2. CML530 REFERENCE"
[[ -s "${CML_GZ}.fai" ]] || samtools faidx "$CML_GZ"
printf 'CML530 sequences: %s | total bp: %s\n' \
  "$(wc -l < "${CML_GZ}.fai")" "$(awk '{s+=$2} END{print s}' "${CML_GZ}.fai")"
printf 'first 12 seq names: %s\n' "$(cut -f1 "${CML_GZ}.fai" | head -12 | tr '\n' ' ')"
printf 'B73 v5 first 12   : %s\n' "$(cut -f1 "${B73V5_FA}.fai" | head -12 | tr '\n' ' ')"

# ---- 3. STEP 1: tags -> CML530 ----------------------------------------------
rule "3. MAP TAGS -> CML530"
# indexes are built once by agent/qtl_index_genomes.sh (~35 min); fail loudly
# rather than silently starting a 35 min build inside the pipeline.
if [[ ! -s "${CML_IDX}.mbw" ]]; then
  printf 'ERROR: %s.mbw missing. Run agent/qtl_index_genomes.sh first.\n' "$CML_IDX"
  exit 1
fi
minibwa map -t "$THREADS" "$CML_IDX" "$TAGS" 2>"${OUT}/map1.stderr" \
  | samtools sort -@2 -o "${OUT}/tags_vs_${ACCESSION}.bam"
samtools index "${OUT}/tags_vs_${ACCESSION}.bam"
printf 'alignment summary (tags -> %s):\n' "$ACCESSION"
samtools flagstat "${OUT}/tags_vs_${ACCESSION}.bam"
printf '\nMAPQ distribution:\n'
samtools view "${OUT}/tags_vs_${ACCESSION}.bam" | awk '{print $5}' \
  | sort -n | uniq -c | awk '{printf "  MAPQ %-4s %s\n", $2, $1}' | tail -15

# ---- 4. STEP 2: extract aligned span, map -> B73 v5 -------------------------
rule "4. EXTRACT CML530 SPAN -> MAP TO B73 v5"
# aligned span only (no flank). BED of the CML530 reference footprint per tag.
samtools view -F 0x904 -q 20 "${OUT}/tags_vs_${ACCESSION}.bam" \
  | awk 'BEGIN{OFS="\t"}{
      n=0; c=$6; while (match(c, /^[0-9]+[MIDNSHP=X]/)) {
        L=substr(c,RSTART,RLENGTH-1); O=substr(c,RSTART+RLENGTH-1,1);
        if (O ~ /[MDN=X]/) n+=L; c=substr(c,RSTART+RLENGTH) }
      print $3, $4-1, $4-1+n, $1
    }' > "${OUT}/${ACCESSION}_spans.bed"
printf 'spans extracted: %s\n' "$(wc -l < "${OUT}/${ACCESSION}_spans.bed")"
# samtools faidx -r, not bedtools getfasta: samtools reads BGZF directly, bedtools
# would need the 2.3 GB decompressed copy. Region list is 1-based inclusive, so
# convert from the 0-based BED start. Read name = the tag's AlleleID, preserved so
# a span can still be traced back to its marker for diagnostics (the GENOTYPE join
# is by v5 coordinate, not by name).
awk 'BEGIN{OFS=""}{print $1,":",$2+1,"-",$3}' "${OUT}/${ACCESSION}_spans.bed" \
  > "${OUT}/${ACCESSION}_spans.regions"
cut -f4 "${OUT}/${ACCESSION}_spans.bed" > "${OUT}/${ACCESSION}_spans.names"
# CRITICAL: rename each record to its unique AlleleID. samtools faidx -r names
# records by REGION, and several tags map to the SAME CML530 span, producing
# duplicate FASTA names. minibwa treats consecutive same-named records as a mate
# pair and hard-aborts on a third:
#   Assertion failed: (i - j <= 2), function worker_pipeline, map-main.c:144
# faidx emits one record per region line, in order, so a positional rename is
# safe. The original coordinate is kept as a trailing comment (aligners use only
# the first whitespace-delimited token as the name).
samtools faidx -r "${OUT}/${ACCESSION}_spans.regions" "$CML_GZ" \
  | awk -v nf="${OUT}/${ACCESSION}_spans.names" '
      BEGIN { while ((getline l < nf) > 0) nm[++k] = l }
      /^>/  { i++; print ">" nm[i] " " substr($0, 2); next }
              { print }' \
  > "${OUT}/${ACCESSION}_tagspans.fa"
nspan=$(grep -c '^>' "${OUT}/${ACCESSION}_tagspans.fa")
ndup=$(grep '^>' "${OUT}/${ACCESSION}_tagspans.fa" | awk '{print $1}' | sort | uniq -d | wc -l | tr -d ' ')
printf 'span sequences: %s | duplicate names after rename: %s (must be 0)\n' "$nspan" "$ndup"
[[ "$ndup" -eq 0 ]] || { printf 'ERROR: duplicate span names remain\n'; exit 1; }
printf 'tags dropped before extraction (secondary/supplementary or MAPQ<20): %s\n' \
  "$(( $(samtools view -c "${OUT}/tags_vs_${ACCESSION}.bam") - nspan ))"

if [[ ! -s "${B73_IDX}.mbw" ]]; then
  printf 'ERROR: %s.mbw missing. Run agent/qtl_index_genomes.sh first.\n' "$B73_IDX"
  exit 1
fi
minibwa map -t "$THREADS" "$B73_IDX" "${OUT}/${ACCESSION}_tagspans.fa" \
  2>"${OUT}/map2.stderr" \
  | samtools sort -@2 -o "${OUT}/${ACCESSION}_vs_B73v5.bam"
samtools index "${OUT}/${ACCESSION}_vs_B73v5.bam"
printf 'alignment summary (%s spans -> B73 v5):\n' "$ACCESSION"
samtools flagstat "${OUT}/${ACCESSION}_vs_B73v5.bam"

# ---- 5. STEP 3: genotype positionally at the marker sites -------------------
rule "5. PILEUP AT MARKER v5 POSITIONS"
bcftools mpileup -f "$B73V5_FA" -T "$BED" -a AD,DP -Ou \
    "${OUT}/${ACCESSION}_vs_B73v5.bam" \
  | bcftools call -m -Ov -o "${OUT}/${ACCESSION}_at_markers.vcf"
printf 'VCF records: %s\n' "$(bcftools view -H "${OUT}/${ACCESSION}_at_markers.vcf" | wc -l)"
# bcftools query confirmed present in bcftools 1.23.1
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t[%GT]\t[%DP]\n' \
  "${OUT}/${ACCESSION}_at_markers.vcf" > "${OUT}/${ACCESSION}_calls.tsv"
printf 'wrote %s (%s rows)\n' "${OUT}/${ACCESSION}_calls.tsv" \
  "$(wc -l < "${OUT}/${ACCESSION}_calls.tsv")"

# ---- 6. classify each marker against the roster -----------------------------
rule "6. CLASSIFY: ref / alt / no-call / no-hit"
Rscript -e '
  suppressMessages(library(data.table))
  a <- commandArgs(TRUE)
  roster <- fread(a[1], colClasses = list(character = "chr_v5"))
  calls  <- fread(a[2], header = FALSE, col.names =
                  c("chrom", "pos", "ref", "alt", "qual", "gt", "dp"))
  roster[, chrom := paste0("chr", chr_v5)]
  # join on the v5 COORDINATE -- the whole point of routing through v5
  m <- merge(roster, calls, by.x = c("chrom", "pos_v5"), by.y = c("chrom", "pos"),
             all.x = TRUE)
  m[, cml_class := fifelse(is.na(gt), "no-hit",
                    fifelse(alt == "." | alt == "", "ref",
                     fifelse(gt %in% c("1/1", "1|1"), "alt",
                      fifelse(gt %in% c("0/1", "0|1", "1/0"), "het-call", "other"))))]
  cat("\n-- all markers --\n");        print(m[, .N, by = cml_class][order(-N)])
  cat("\n-- framework markers --\n");  print(m[framework == TRUE, .N, by = cml_class][order(-N)])
  cat("\n-- by chromosome (framework) --\n")
  print(dcast(m[framework == TRUE], chr_v5 ~ cml_class, value.var = "marker",
              fun.aggregate = length)[order(as.integer(chr_v5))])
  cat("\nNOTE: het-call in an inbred line is a red flag -- it means the span\n")
  cat("aligned to a paralog or a collapsed repeat, not that CML530 is truly\n")
  cat("heterozygous. Treat het-call as unresolved, NOT as a genotype.\n")
  out <- m[, .(marker, CloneID, chr_v5, pos_v5, chr_v4, pos_v4, SNP, SnpPosition,
               b73_ref = ref, cml_alt = alt, gt, dp, qual, cml_class,
               qc_pass, framework)]
  setorder(out, chr_v5, pos_v5)
  fwrite(out, a[3], sep = "\t")
  cat("\nwrote", a[3], "rows:", nrow(out), "\n")
' "$ROSTER" "${OUT}/${ACCESSION}_calls.tsv" "${DERIVED}/${ACCESSION}_marker_alleles.tsv"

rule "DONE"
printf 'outputs in %s\n' "$OUT"
ls -la "$OUT"
printf '\nfinal table: %s\n' "${DERIVED}/${ACCESSION}_marker_alleles.tsv"
printf '\nINTERPRETATION\n'
printf '  1. CML530 is a FEMALE PARENT of this population (README.md:44-46), so\n'
printf '     these are actual parental alleles. ABH parent A := the CML530\n'
printf '     allele; parent B := the complement, which is sound because the\n'
printf '     markers were filtered to those segregating 1:2:1 (i.e. the parents\n'
printf '     are homozygous for different alleles).\n'
printf '  2. "ref" means CML530 matches B73 -- DArT defined the ref allele by its\n'
printf '     AGPv4 alignment. It says nothing about the ENxx tester. Polarity is\n'
printf '     PER MARKER: CML530 carries the B73 allele at some sites and the\n'
printf '     alternate at others, with no global pattern.\n'
printf '  3. STILL ASSUMED, and to confirm with Stelio: that this F2 came from a\n'
printf '     CML530 x ENxx cross, and WHICH tester. paper/maize_pue_manuscript.md:188\n'
printf '     lists CML530 x EN31 among the most stable high-yielding hybrids by\n'
printf '     WAASY -- a lead, not a confirmation (stability != P-contrast).\n'
printf '  4. het-call in an inbred parent is NOT heterozygosity; it means the tag\n'
printf '     span hit a paralog or collapsed repeat. Treated as unresolved.\n'
printf '  5. no-hit rate should rise toward pericentromeres (expected); a no-hit\n'
printf '     cluster on a distal arm is a result worth looking at, not noise.\n'
printf '     The functional test already showed 88 of 3828 tags with no CML530\n'
printf '     locus (97.7%% mapped) -- candidate presence/absence variants.\n'
