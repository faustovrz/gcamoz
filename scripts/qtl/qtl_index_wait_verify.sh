#!/usr/bin/env bash
# qtl_index_wait_verify.sh -- wait for the in-flight minibwa index to finish, then
# verify BOTH indexes are complete and usable.
#
# WHY THIS EXISTS: the background shell wrapping agent/qtl_index_genomes.sh was
# torn down mid-run (no completion record). The `minibwa index` process SURVIVED
# and kept writing, but the wrapper's own logging and its "DONE" summary never
# execute -- so qtl_index_genomes.log is permanently truncated after
# "INDEX B73 v5" and cannot be used to confirm success.
#
# Therefore completeness is established from the ARTEFACTS, not the log:
#   1. no `minibwa index` process remains
#   2. both .l2b and .mbw exist for each genome
#   3. .l2b size == ceil(genome_bp / 4) +/- slack  (2-bit packing; this is the
#      signature of a COMPACTED, finished .l2b -- a partial one is ~2x larger,
#      which is why B73v5.l2b read 1.09 GB mid-run and 545 MB when done)
#   4. FUNCTIONAL test: minibwa map actually aligns the real tag FASTA against
#      each index and returns a plausible number of alignments. A truncated .mbw
#      passes a size check but fails here.
#
# Usage (background):
#   bash agent/qtl_index_wait_verify.sh > agent/qtl_index_wait_verify.log 2>&1

set -uo pipefail

ZEAL="/Users/fvrodriguez/repos/zealhmm"
IDXDIR="${ZEAL}/data/ref/genomes/minibwa_index"
CML_GZ="${ZEAL}/data/ref/genomes/Zm-CML530-REFERENCE-HiLo-1.0/Zm-CML530-REFERENCE-HiLo-1.0.fa.gz"
B73V5_FA="/Users/fvrodriguez/repos/zealtiger/data/brbseq/ref/Zm-B73-REFERENCE-NAM-5.0.fa"
TAGS="/Users/fvrodriguez/repos/gcamoz/data/qtl/derived/tags_for_cml530.fa"

rule() { printf '\n%s\n%s\n%s\n' "======================================================================" "$1" "======================================================================"; }

rule "0. WAIT FOR minibwa index TO EXIT"
printf 'start waiting: %s\n' "$(date +"%Y-%m-%d %H:%M:%S %Z")"
n=0
until ! pgrep -f "minibwa index" >/dev/null 2>&1; do
  sz=$(stat -f%z "${IDXDIR}/B73v5.mbw" 2>/dev/null || echo 0)
  printf '  [%s] B73v5.mbw = %.2f GB\n' "$(date +%H:%M:%S)" "$(echo "$sz" | awk '{print $1/1073741824}')"
  n=$((n+1)); [[ $n -gt 240 ]] && { printf 'TIMEOUT after ~40 min of waiting\n'; break; }
  sleep 30
done
printf 'no minibwa index process remaining: %s\n' "$(date +"%Y-%m-%d %H:%M:%S %Z")"
sleep 3   # let the final write flush

rule "1. ARTEFACT CHECK"
fail=0
check_one() {
  local name="$1" pfx="$2" fa="$3" faidx="$4"
  local bp exp l2b mbw
  bp=$(awk '{s+=$2} END{print s}' "$faidx")
  exp=$(echo "$bp" | awk '{printf "%d", $1/4}')
  l2b=$(stat -f%z "${pfx}.l2b" 2>/dev/null || echo 0)
  mbw=$(stat -f%z "${pfx}.mbw" 2>/dev/null || echo 0)
  printf '\n%s\n' "$name"
  printf '  genome bp        : %s\n' "$bp"
  printf '  .l2b expected    : ~%s bytes (bp/4, 2-bit packed)\n' "$exp"
  printf '  .l2b actual      : %s\n' "$l2b"
  printf '  .mbw actual      : %s (%.2f GB)\n' "$mbw" "$(echo "$mbw" | awk '{print $1/1073741824}')"
  # .l2b within 5% of bp/4 => compacted/finished
  local ok_l2b
  ok_l2b=$(echo "$l2b $exp" | awk '{print ($1 > $2*0.95 && $1 < $2*1.05) ? 1 : 0}')
  if [[ "$ok_l2b" -eq 1 ]]; then printf '  .l2b             : OK (compacted)\n'
  else printf '  .l2b             : SUSPECT -- not the compacted size\n'; fail=1; fi
  if [[ "$mbw" -gt 1000000000 ]]; then printf '  .mbw             : present\n'
  else printf '  .mbw             : SUSPECT -- too small\n'; fail=1; fi
}
[[ -s "${CML_GZ}.fai" ]] || samtools faidx "$CML_GZ"
check_one "CML530" "${IDXDIR}/CML530" "$CML_GZ"   "${CML_GZ}.fai"
check_one "B73 v5" "${IDXDIR}/B73v5"  "$B73V5_FA" "${B73V5_FA}.fai"

rule "2. FUNCTIONAL TEST (the check a size test cannot make)"
if [[ ! -s "$TAGS" ]]; then
  printf 'tag FASTA missing: %s -- skipping functional test\n' "$TAGS"; fail=1
else
  ntags=$(grep -c '^>' "$TAGS")
  printf 'tags: %s\n' "$ntags"
  for g in CML530 B73v5; do
    printf '\n-- map tags -> %s --\n' "$g"
    if out=$(minibwa map -t 4 "${IDXDIR}/${g}" "$TAGS" 2>/dev/null); then
      tot=$(echo "$out" | grep -vc '^@')
      mapped=$(echo "$out" | awk '!/^@/ && $3!="*"' | wc -l | tr -d ' ')
      q30=$(echo "$out" | awk '!/^@/ && $5>=30' | wc -l | tr -d ' ')
      printf '  records %s | mapped %s | MAPQ>=30 %s\n' "$tot" "$mapped" "$q30"
      if [[ "$mapped" -lt $((ntags / 2)) ]]; then
        printf '  SUSPECT -- fewer than half the tags aligned\n'; fail=1
      else printf '  OK\n'; fi
    else
      printf '  FAILED to run minibwa map against %s\n' "$g"; fail=1
    fi
  done
fi

rule "VERDICT"
printf 'total index dir: %s\n' "$(du -sh "$IDXDIR" | cut -f1)"
printf 'free disk     : %s\n' "$(df -h "$IDXDIR" | awk 'NR==2{print $4}')"
if [[ "$fail" -eq 0 ]]; then
  printf '\nBOTH INDEXES COMPLETE AND FUNCTIONAL. Safe to run\n'
  printf 'agent/qtl_cml530_genotypes.sh\n'
else
  printf '\nPROBLEM DETECTED -- do NOT run the pipeline. Re-run\n'
  printf 'agent/qtl_index_genomes.sh (it skips whichever index is already valid).\n'
fi
printf 'finished: %s\n' "$(date +"%Y-%m-%d %H:%M:%S %Z")"
exit "$fail"
