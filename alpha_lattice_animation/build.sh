#!/usr/bin/env bash
# Render all 14 scenes (S01-S14) and concatenate them into one narrated film.
# Usage:  ./build.sh            (low quality, fast preview)
#         ./build.sh -qh        (1080p60 final)
set -euo pipefail
cd "$(dirname "$0")"

QUAL="${1:--ql}"                        # -ql (480p15) default; pass -qh for 1080p60
MANIM=.venv/bin/manim
SCENES=(S01_Optimize S02_Problem S03_GeneratingArray S04_CyclicDevelopment S05_RowOffsets \
        S06_Blocks S07_CoocPart1 S08_Randomization S09_Alpha01Design S10_Alpha01Build S11_SingleCoOccur \
        S12_CoocPart2 S13_RealTrial S14_CoocReal)

# 1. render every scene (gTTS narration is cached after first run)
"$MANIM" "$QUAL" alpha_lattice.py "${SCENES[@]}"

# 2. locate the per-scene mp4s (quality subdir name depends on -ql/-qh)
QDIR=$(ls -d media/videos/alpha_lattice/*/ | head -1)
mkdir -p output
LIST=output/_concat.txt
: > "$LIST"
for s in "${SCENES[@]}"; do
  echo "file '$(cd "$QDIR" && pwd)/$s.mp4'" >> "$LIST"
done

# 3. concat via the demuxer (correct per-scene durations), re-encoding the audio into one
#    continuous AAC stream (removes the AAC join pops that `-c copy` leaves). The
#    `aresample=async=1` filter is ESSENTIAL: each scene's audio is shorter than its video
#    (held final frame + trailing pause), and without it the re-encoder collapses those
#    inter-scene silences, drifting the audio ahead of the picture and clipping later
#    narration. async=1 fills the gaps with silence so audio and video stay locked.
ffmpeg -y -f concat -safe 0 -i "$LIST" \
  -c:v copy -af aresample=async=1 -c:a aac -b:a 192k output/alpha_lattice_full.mp4
rm -f "$LIST"
echo "Done -> output/alpha_lattice_full.mp4"
