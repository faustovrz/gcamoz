#!/usr/bin/env bash
# install_minibwa_from_source.sh -- build minibwa WITH the GPL low-memory BWT.
#
# WHY NOT HOMEBREW: `brew install minibwa` (0.6) is built with gpl=0, so the
# low-memory BWT construction path is absent:
#     $ minibwa index -l -t 8 chr10.fa c10l
#     ERROR: option -l not compiled as it depends on GPL'd code
# The -l flag still appears in --help, so the help text is NOT evidence that it
# works. Verified failing 2026-07-30.
#
# We need -l: minibwa index uses 18N RAM (N = genome size in bp). Maize is
# ~2.3 Gb -> ~42 GB, and this machine has 24 GB. Default `make` includes the
# GPL path (the Makefile only drops it when gpl=0 is passed), so a source build
# is the fix.
#
# Installs to ~/bio/bin, which precedes /opt/homebrew/bin on PATH and therefore
# shadows the Homebrew binary.
#
# Usage: bash agent/install_minibwa_from_source.sh > agent/install_minibwa_from_source.log 2>&1

set -euo pipefail

SRC="${HOME}/bio/src"
BIN="${HOME}/bio/bin"
REPO_URL="https://github.com/lh3/minibwa"

rule() { printf '\n%s\n%s\n%s\n' "======================================================================" "$1" "======================================================================"; }

rule "0. ENVIRONMENT"
printf 'arch     : %s\n' "$(uname -m)"
printf 'cc       : %s\n' "$(cc --version 2>&1 | head -1)"
printf 'target   : %s\n' "$BIN"
printf 'brew ver : %s\n' "$(brew list --versions minibwa 2>/dev/null || echo 'not installed')"

rule "1. CLONE"
mkdir -p "$SRC" "$BIN"
if [[ -d "${SRC}/minibwa/.git" ]]; then
  git -C "${SRC}/minibwa" fetch --tags origin
  git -C "${SRC}/minibwa" checkout master
  git -C "${SRC}/minibwa" pull --ff-only
else
  git clone "$REPO_URL" "${SRC}/minibwa"
fi
cd "${SRC}/minibwa"
printf 'commit: %s\n' "$(git rev-parse --short HEAD)"
printf 'date  : %s\n' "$(git log -1 --format=%cI)"

rule "2. BUILD (default flags -> GPL low-memory BWT INCLUDED)"
make clean >/dev/null 2>&1 || true
# no gpl=0, no mimalloc=0: we want the full build. omp is auto-detected.
make -j"$(sysctl -n hw.ncpu)" 2>&1 | tail -15
printf '\nbuilt binary:\n'
ls -la ./minibwa

rule "3. VERIFY -l IS PRESENT"
printf 'version: %s\n' "$(./minibwa version)"
TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT
printf '>t\n%s\n' "$(head -c 200000 /dev/urandom | base64 | tr -dc 'ACGT' | head -c 100000)" > "${TMP}/t.fa"
if ./minibwa index -l -t 4 "${TMP}/t.fa" "${TMP}/t" >/dev/null 2>"${TMP}/err"; then
  printf 'RESULT: -l WORKS (low-memory BWT available)\n'
else
  printf 'RESULT: -l STILL FAILS:\n'; cat "${TMP}/err"; exit 1
fi

rule "4. INSTALL"
install -m 0755 ./minibwa "${BIN}/minibwa"
printf 'installed: %s\n' "${BIN}/minibwa"
hash -r 2>/dev/null || true
printf 'which minibwa -> %s\n' "$(command -v minibwa)"
printf 'version      -> %s\n' "$(minibwa version)"
printf '\nNOTE: %s precedes /opt/homebrew/bin on PATH, so this shadows the\n' "$BIN"
printf 'Homebrew 0.6 bottle. Consider `brew uninstall minibwa` to avoid ambiguity.\n'
