#!/usr/bin/env bash
# qtl_minibwa_index_bench.sh -- measure minibwa index cost before committing to
# indexing two ~2.2 Gb maize genomes on a 24 GB box.
#
# Indexes B73 v5 chromosome 10 (~152 Mb, the smallest chromosome) with the same
# flags the real pipeline uses (-l low-memory BWT, -t THREADS), times it, and
# reports on-disk index sizes. libsais SA construction is near-linear in genome
# size, so scaling chr10 -> whole genome by bp ratio gives a usable estimate for
# both time and index size.
#
# Also reports the uncompressed size of the CML530 assembly (gzip -l, instant)
# and times the gunzip step, which is on the pipeline's critical path.
#
# Usage: bash agent/qtl_minibwa_index_bench.sh > agent/qtl_minibwa_index_bench.log 2>&1

set -euo pipefail

ZEAL="/Users/fvrodriguez/repos/zealhmm"
B73V5_FA="/Users/fvrodriguez/repos/zealtiger/data/brbseq/ref/Zm-B73-REFERENCE-NAM-5.0.fa"
CML_GZ="${ZEAL}/data/ref/genomes/Zm-CML530-REFERENCE-HiLo-1.0/Zm-CML530-REFERENCE-HiLo-1.0.fa.gz"
BENCH="$(mktemp -d)"
THREADS="${THREADS:-8}"
trap 'rm -rf "$BENCH"' EXIT

rule() { printf '\n%s\n%s\n%s\n' "======================================================================" "$1" "======================================================================"; }

rule "0. GENOME SIZES"
B73_BP=$(awk '{s+=$2} END{print s}' "${B73V5_FA}.fai")
printf 'B73 v5      : %s bp  (%.2f Gb) across %s sequences\n' \
  "$B73_BP" "$(echo "$B73_BP" | awk '{print $1/1e9}')" "$(wc -l < "${B73V5_FA}.fai")"
# CML530 ships as BGZF (gzip-compatible, block-compressed), so gzip -l reports 0.
# samtools can faidx it in place -- no 2.2 GB decompression needed anywhere.
[[ -s "${CML_GZ}.fai" ]] || samtools faidx "$CML_GZ"
CML_BP=$(awk '{s+=$2} END{print s}' "${CML_GZ}.fai")
printf 'CML530      : %s bp  (%.2f Gb) across %s sequences [BGZF, faidx in place]\n' \
  "$CML_BP" "$(echo "$CML_BP" | awk '{print $1/1e9}')" "$(wc -l < "${CML_GZ}.fai")"
printf 'CML530 seq names (first 12): %s\n' "$(cut -f1 "${CML_GZ}.fai" | head -12 | tr '\n' ' ')"
printf 'threads     : %s | RAM: %s GB\n' \
  "$THREADS" "$(( $(sysctl -n hw.memsize) / 1024 / 1024 / 1024 ))"

rule "1. BENCHMARK: index B73 v5 chr10 with production flags"
# NB: this B73 v5 FASTA names chromosomes chr1..chr10, NOT 1..10 (see notes).
samtools faidx "$B73V5_FA" chr10 > "${BENCH}/chr10.fa"
CHR10_BP=$(grep -v '^>' "${BENCH}/chr10.fa" | tr -d '\n' | wc -c | tr -d ' ')
printf 'chr10: %s bp (%.1f Mb)\n' "$CHR10_BP" "$(echo "$CHR10_BP" | awk '{print $1/1e6}')"
printf 'running: minibwa index -l -t %s chr10.fa\n\n' "$THREADS"

T0=$(date +%s)
/usr/bin/time -l minibwa index -l -t "$THREADS" "${BENCH}/chr10.fa" "${BENCH}/chr10" \
  2>&1 | grep -Ei "real|maximum resident|CMD|Real time|Peak RSS" || true
T1=$(date +%s)
ELAPSED=$(( T1 - T0 ))
printf '\nwall clock: %s s\n' "$ELAPSED"

printf '\nindex files produced:\n'
ls -la "${BENCH}"/chr10.* | awk '{printf "  %-12s %s\n", $5, $NF}'
IDX_BYTES=$(du -ck "${BENCH}"/chr10.l2b "${BENCH}"/chr10.mbw 2>/dev/null | tail -1 | awk '{print $1*1024}')
printf 'total index bytes: %s (%.2f bytes per bp)\n' \
  "$IDX_BYTES" "$(echo "$IDX_BYTES $CHR10_BP" | awk '{print $1/$2}')"

rule "2. EXTRAPOLATION TO FULL GENOMES"
awk -v e="$ELAPSED" -v c="$CHR10_BP" -v ib="$IDX_BYTES" \
    -v b="$B73_BP" -v m="$CML_BP" 'BEGIN{
  bps = c / e;  bpb = ib / c;
  printf "measured throughput : %.1f Mb/s   index density: %.2f bytes/bp\n\n", bps/1e6, bpb;
  printf "%-12s %12s %10s %14s\n", "genome", "bp", "index_time", "index_size";
  printf "%-12s %12d %8.1f min %10.2f GB\n", "B73 v5",  b, (b/bps)/60, (b*bpb)/1e9;
  printf "%-12s %12d %8.1f min %10.2f GB\n", "CML530",  m, (m/bps)/60, (m*bpb)/1e9;
  printf "%-12s %12d %8.1f min %10.2f GB\n", "TOTAL", b+m, ((b+m)/bps)/60, ((b+m)*bpb)/1e9;
}'

rule "3. CML530 DECOMPRESSION (critical path, one-off)"
printf 'timing gunzip of 200 MB sample to get throughput...\n'
T0=$(date +%s%N)
gunzip -c "$CML_GZ" 2>/dev/null | head -c 200000000 > /dev/null || true
T1=$(date +%s%N)
GZ_MS=$(( (T1 - T0) / 1000000 ))
awk -v ms="$GZ_MS" -v m="$CML_BP" 'BEGIN{
  mbs = 200.0 / (ms/1000.0);
  printf "decompress throughput: %.0f MB/s -> full %.2f Gb takes ~%.1f min, %.2f GB on disk\n",
         mbs, m/1e9, (m/1e6/mbs)/60, m/1e9;
}'

rule "NOTES"
cat <<'EOF'
- Mapping itself is negligible: 3,828 tags of 36-69 bp against a prebuilt index
  is seconds, not minutes. Indexing is essentially the entire cost.
- Indexes are built ONCE and reused. Adding another accession later (the ENxx
  tester, or another CML female) costs one more index each; B73 v5 is permanent.
  NB: CML530 is itself a parent of this population, not a stand-in.
- -l trades speed for memory. Without it minibwa index wants 18N ~= 40 GB,
  which this 24 GB box does not have, so -l is not optional here.
- Extrapolation assumes near-linear libsais SA construction. Whole-genome runs
  usually come in somewhat WORSE than a single-chromosome extrapolation because
  of cache pressure and repeat content, so treat these as a floor.
EOF
