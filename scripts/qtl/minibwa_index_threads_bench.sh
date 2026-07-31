#!/usr/bin/env bash
# qtl_minibwa_index_threads_bench.sh -- does minibwa index scale with threads?
#
# CORRECTS AN ERROR in agent/qtl_minibwa_index_bench.sh. That benchmark ran with
# -l (the low-memory GPL BWT) and observed 68.97s CPU vs 69.11s wall on -t 8,
# i.e. no thread scaling. I then used that to argue an HPC cluster offers no
# benefit. That was circular: -l is required ONLY because this laptop has 24 GB
# and the default path wants 18N (~42 GB for maize). A big-memory node would not
# use -l, and the default path uses libsais with OpenMP, which DOES parallelize.
#
# So the honest question is: how fast is the DEFAULT (non -l) path, and does it
# scale with -t? chr10 is 152 Mb -> 18N ~= 2.7 GB, which fits in 24 GB, so the
# parallel path can be measured here directly and extrapolated to a large-memory
# node.
#
# Measures, all on B73 v5 chr10:
#   (a) -l        -t 8   (the laptop-constrained path, for reference)
#   (b) default   -t 1   (serial baseline of the fast path)
#   (c) default   -t 8   (does it scale?)
#
# Usage: bash agent/qtl_minibwa_index_threads_bench.sh > agent/qtl_minibwa_index_threads_bench.log 2>&1

set -euo pipefail

B73V5_FA="/Users/fvrodriguez/repos/zealtiger/data/brbseq/ref/Zm-B73-REFERENCE-NAM-5.0.fa"
BENCH="$(mktemp -d)"
trap 'rm -rf "$BENCH"' EXIT

rule() { printf '\n%s\n%s\n%s\n' "======================================================================" "$1" "======================================================================"; }

rule "0. SETUP"
printf 'minibwa : %s (%s)\n' "$(command -v minibwa)" "$(minibwa version)"
printf 'cores   : %s | RAM: %s GB\n' "$(sysctl -n hw.ncpu)" "$(( $(sysctl -n hw.memsize) / 1024 / 1024 / 1024 ))"
samtools faidx "$B73V5_FA" chr10 > "${BENCH}/chr10.fa"
BP=$(grep -v '^>' "${BENCH}/chr10.fa" | tr -d '\n' | wc -c | tr -d ' ')
printf 'chr10   : %s bp (%.1f Mb)  -> default path needs 18N ~= %.1f GB\n' \
  "$BP" "$(echo "$BP" | awk '{print $1/1e6}')" "$(echo "$BP" | awk '{print 18*$1/1e9}')"

# run one configuration; echo "label elapsed_s peak_gb"
run_one() {
  local label="$1"; shift
  local pfx="${BENCH}/idx_$(echo "$label" | tr -dc 'a-z0-9')"
  rm -f "${pfx}".* 2>/dev/null || true
  local t0 t1 out
  t0=$(date +%s%N)
  out=$(/usr/bin/time -l minibwa index "$@" "${BENCH}/chr10.fa" "$pfx" 2>&1)
  t1=$(date +%s%N)
  local ms=$(( (t1 - t0) / 1000000 ))
  local rss cpu
  rss=$(echo "$out" | awk '/maximum resident set size/{print $1}')
  cpu=$(echo "$out" | grep -oE 'CPU: [0-9.]+' | awk '{print $2}')
  printf '%-22s wall %7.1f s   CPU %8s s   peakRSS %6.2f GB\n' \
    "$label" "$(echo "$ms" | awk '{print $1/1000}')" "${cpu:-NA}" \
    "$(echo "${rss:-0}" | awk '{print $1/1073741824}')"
  echo "$ms" > "${BENCH}/last_ms_${label//[^a-z0-9]/}"
}

rule "1. MEASUREMENTS (B73 v5 chr10)"
run_one "lowmem -l -t8"   -l -t 8
run_one "default -t1"        -t 1
run_one "default -t8"        -t 8

rule "2. INTERPRETATION"
L=$(cat "${BENCH}/last_ms_lowmemlt8"); D1=$(cat "${BENCH}/last_ms_defaultt1"); D8=$(cat "${BENCH}/last_ms_defaultt8")
awk -v l="$L" -v d1="$D1" -v d8="$D8" -v bp="$BP" 'BEGIN{
  printf "thread speedup of default path (t1/t8) : %.2fx\n", d1/d8;
  printf "default -t8 vs low-memory -l -t8       : %.2fx faster\n", l/d8;
  printf "\nextrapolated to a 2.18 Gb genome (linear in bp):\n";
  s = 2182075994.0/bp;
  printf "  %-24s %6.1f min\n", "laptop  (-l -t8)",    (l*s/1000)/60;
  printf "  %-24s %6.1f min\n", "big-mem (default t1)",(d1*s/1000)/60;
  printf "  %-24s %6.1f min\n", "big-mem (default t8)",(d8*s/1000)/60;
  printf "\nNOTE the default path needs 18N ~= 39 GB for a 2.18 Gb genome, so the\n";
  printf "big-mem rows are NOT achievable on this 24 GB laptop -- they are what a\n";
  printf "cluster node would deliver. On Hazel more -t would likely go further\n";
  printf "than 8 threads, so treat the t8 row as an upper bound on time.\n";
}'
