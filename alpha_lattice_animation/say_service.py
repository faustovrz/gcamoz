"""A manim-voiceover SpeechService backed by the macOS `say` command.

Local, offline, no API key, and far more natural than gTTS. Each line is synthesized
with `say -v <voice>` to AIFF, then converted to MP3 (the format the rest of the
manim-voiceover pipeline expects). Caching mirrors the built-in GTTSService.
"""

from __future__ import annotations

import re
import shutil
import subprocess
from pathlib import Path

from manim import logger

from manim_voiceover.helper import remove_bookmarks
from manim_voiceover.services.base import (
    PathLike,
    SpeechService,
    initialize_speech_service,
    path_to_string,
)


class SayService(SpeechService):
    """Text-to-speech via the macOS `say` binary (NSSpeechSynthesizer voices)."""

    def __init__(self, voice: str = "Evan (Enhanced)", rate: int | None = 155,
                 sentence_pause: int = 550, comma_pause: int = 180,
                 dash_pause: int = 320, trailing_pause: int = 750,
                 letter_pause: int = 250, **kwargs: object) -> None:
        """
        Args:
            voice: a `say` voice name. List them with `say -v '?'`. "(Enhanced)" /
                "(Premium)" voices must first be downloaded in System Settings.
            rate:  words per minute; None = voice default (~175). Lower = slower.
            sentence_pause: milliseconds of silence inserted after . ? !
            comma_pause:    milliseconds of silence inserted after commas.
            dash_pause:     milliseconds of silence for " -- " breaks.
            trailing_pause: milliseconds of silence appended to the END of every clip,
                giving a clear gap between sections and between scenes.
            letter_pause:   milliseconds of silence BEFORE a standalone capital-letter word
                (e.g. the factor "E"), so `say` does not glue it onto the previous word
                ("factor E" -> "factory").
        """
        if shutil.which("say") is None:
            raise RuntimeError("macOS `say` not found — SayService only runs on macOS.")
        if shutil.which("ffmpeg") is None:
            raise RuntimeError("ffmpeg not found — needed to convert `say` output to mp3.")
        initialize_speech_service(self, kwargs)
        self.voice = voice
        self.rate = rate
        self.sentence_pause = sentence_pause
        self.comma_pause = comma_pause
        self.dash_pause = dash_pause
        self.trailing_pause = trailing_pause
        self.letter_pause = letter_pause

    def _with_pauses(self, text: str) -> str:
        """Insert `say` silence commands so punctuation gets real, audible pauses,
        and append a trailing silence so each section/scene ends with a clear gap."""
        t = re.sub(r"\s+--\s+", f" [[slnc {self.dash_pause}]] ", text)
        if self.letter_pause:
            # brief gap before a standalone capital letter (e.g. "factor E") so `say`
            # pronounces it as a separate word instead of slurring it ("factory").
            t = re.sub(r"(?<=\w )([A-Z])\b", rf"[[slnc {self.letter_pause}]] \1", t)
        if self.comma_pause:
            t = re.sub(r",(\s+)", rf", [[slnc {self.comma_pause}]]\1", t)
        if self.sentence_pause:
            t = re.sub(r"([.?!])(\s+)", rf"\1 [[slnc {self.sentence_pause}]]\2", t)
        if self.trailing_pause:
            t = f"{t} [[slnc {self.trailing_pause}]]"
        return t

    def generate_from_text(
        self,
        text: str,
        cache_dir: PathLike | None = None,
        path: PathLike | None = None,
        **kwargs: object,
    ) -> dict:
        if cache_dir is None:
            cache_dir = self.cache_dir

        input_text = remove_bookmarks(text)
        input_data = {
            "input_text": input_text,
            "service": "macos_say",
            "voice": self.voice,
            "rate": self.rate,
            "sentence_pause": self.sentence_pause,
            "comma_pause": self.comma_pause,
            "dash_pause": self.dash_pause,
            "trailing_pause": self.trailing_pause,
            "letter_pause": self.letter_pause,
        }

        cached_result = self.get_cached_result(input_data, cache_dir)
        if cached_result is not None:
            return cached_result

        if path is None:
            audio_path = self.get_audio_basename(input_data) + ".mp3"
        else:
            audio_path = path_to_string(path)

        mp3_path = Path(cache_dir) / audio_path
        aiff_path = mp3_path.with_suffix(".aiff")

        say_cmd = ["say", "-v", self.voice]
        if self.rate is not None:
            say_cmd += ["-r", str(self.rate)]
        say_cmd += ["-o", str(aiff_path), self._with_pauses(input_text)]

        try:
            subprocess.run(say_cmd, check=True, capture_output=True)
            subprocess.run(
                ["ffmpeg", "-y", "-i", str(aiff_path),
                 "-codec:a", "libmp3lame", "-qscale:a", "4", str(mp3_path)],
                check=True, capture_output=True,
            )
        except subprocess.CalledProcessError as e:
            logger.error(e.stderr.decode("utf-8", "ignore") if e.stderr else str(e))
            raise
        finally:
            aiff_path.unlink(missing_ok=True)

        return {
            "input_text": text,
            "input_data": input_data,
            "original_audio": audio_path,
        }
