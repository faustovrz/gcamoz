#!/usr/bin/env bash
# qtl_index_genomes.sh -- build the two minibwa indexes the CML530 genotyping
# pipeline needs. Split out of qtl_cml530_genotypes.sh because this is the only
# expensive step (~35 min) and is done ONCE, then reused.
#
# Uses -l (low-memory GPL BWT). Not optional on this 24 GB laptop: the default
# path wants 18N, i.e. ~39 GB for a 2.18 Gb genome. Measured cost of -l is 8.9x
# (agent/qtl_minibwa_index_threads_bench.log) -- accepted deliberately, see
# that log for the local-vs-HPC numbers.
#
# Block size left at the default 10m. A larger -b would trade RAM for speed and
# could recover part of the 8.9x, but it is UNTESTED here, so this run takes the
# known-good path rather than risking a 35 min failure.
#
# Index prefixes go to a shared dir so nothing is written into the source
# reference directories (B73 v5 belongs to zealtiger, CML530 to zealhmm).
#
# Outputs per genome: <prefix>.l2b (2-bit sequence) + <prefix>.mbw (BWT + sampled SA)
# Expect ~5.2 GB (B73 v5) + ~5.5 GB (CML530) = ~10.7 GB total.
#
# Usage: bash agent/qtl_index_genomes.sh > agent/qtl_index_genomes.log 2>&1

set -euo pipefail

ZEAL="/Users/fvrodriguez/repos/zealhmm"
CML_GZ="${ZEAL}/data/ref/genomes/Zm-CML530-REFERENCE-HiLo-1.0/Zm-CML530-REFERENCE-HiLo-1.0.fa.gz"
B73V5_FA="/Users/fvrodriguez/repos/zealtiger/data/brbseq/ref/Zm-B73-REFERENCE-NAM-5.0.fa"
IDXDIR="${ZEAL}/data/ref/genomes/minibwa_index"
THREADS="${THREADS:-8}"

rule() { printf '\n%s\n%s\n%s\n' "======================================================================" "$1" "======================================================================"; }

rule "0. SETUP"
mkdir -p "$IDXDIR"
printf 'minibwa   : %s (%s)\n' "$(command -v minibwa)" "$(minibwa version)"
printf 'index dir : %s\n' "$IDXDIR"
printf 'threads   : %s | RAM: %s GB\n' "$THREADS" "$(( $(sysctl -n hw.memsize) / 1024 / 1024 / 1024 ))"
printf 'free disk : %s\n' "$(df -h "$IDXDIR" | awk 'NR==2{print $4}')"
printf 'started   : %s\n' "$(date +"%Y-%m-%d %H:%M:%S %Z")"

# faidx CML530 in place -- it ships as BGZF, so no 2.3 GB decompression needed.
[[ -s "${CML_GZ}.fai" ]] || samtools faidx "$CML_GZ"

index_one() {
  local name="$1" fa="$2" pfx="$3"
  rule "INDEX ${name}"
  if [[ -s "${pfx}.mbw" && -s "${pfx}.l2b" ]]; then
    printf 'already present, skipping:\n'; ls -la "${pfx}".{l2b,mbw}; return 0
  fi
  printf 'source : %s\n' "$fa"
  printf 'prefix : %s\n' "$pfx"
  local t0 t1
  t0=$(date +%s)
  /usr/bin/time -l minibwa index -l -t "$THREADS" "$fa" "$pfx" 2>&1 \
    | grep -Ei "Real time|Peak RSS|maximum resident|real|CMD" || true
  t1=$(date +%s)
  printf '\nwall clock: %s min %s s\n' "$(( (t1-t0)/60 ))" "$(( (t1-t0)%60 ))"
  printf 'index files:\n'
  ls -la "${pfx}".{l2b,mbw} | awk '{printf "  %12s  %s\n", $5, $NF}'
  printf 'total: %s\n' "$(du -ch "${pfx}".{l2b,mbw} | tail -1 | cut -f1)"
}

index_one "CML530"  "$CML_GZ"   "${IDXDIR}/CML530"
index_one "B73 v5"  "$B73V5_FA" "${IDXDIR}/B73v5"

rule "DONE"
printf 'finished  : %s\n' "$(date +"%Y-%m-%d %H:%M:%S %Z")"
printf 'index dir contents:\n'
ls -la "$IDXDIR"
printf '\ntotal index size: %s\n' "$(du -sh "$IDXDIR" | cut -f1)"
printf 'free disk after : %s\n' "$(df -h "$IDXDIR" | awk 'NR==2{print $4}')"
