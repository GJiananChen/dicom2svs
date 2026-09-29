# dicom2svs

Single-file tool, `dicom2svs.py`. It converts DICOM whole-slide images into Aperio `.svs`
files. Dependencies are declared inline (PEP 723), so always run it with
`uv run dicom2svs.py …` and never pip-install them.

- To convert slides for the user, follow `.claude/skills/dicom2svs/SKILL.md`: dry run first,
  convert after confirmation, then report the summary.
- Code layout, top to bottom: scanning, one slide per folder (`scan`) → pre-checks
  (`check_source`; metadata problems are `fatal=False` warnings, only unconvertible
  files fail) → tile passthrough and pyramid building (`TileSource`, `build_level`,
  `convert`) → OpenSlide verification (`verify`) → summary table and CSV (`print_summary`).
- The 40x level must stay a byte-for-byte copy of the DICOM JPEG tiles. Only the
  lower levels are re-encoded.
- Must run on Windows, Linux and macOS: no POSIX-only calls (e.g. `os.pread`), open
  files in binary mode, write text as UTF-8.
- After editing, run `uv run tests/run_tests.py` (synthetic slides, no patient data). Never commit real slide data.
