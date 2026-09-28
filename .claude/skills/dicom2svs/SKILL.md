---
name: dicom2svs
description: Convert DICOM whole-slide images (pathology WSI, e.g. Leica/Aperio GT450 or Sectra exports) into pyramidal Aperio .svs files, checking the source is 40x full resolution first. Use when the user wants to convert DICOM slides/folders to SVS, find the full-resolution DICOM in a slide folder, or batch-convert WSI DICOM for QuPath/OpenSlide.
---

# dicom2svs

Converts the largest (full-resolution, 40x) DICOM file of each slide into an Aperio `.svs`.
It copies the 40x JPEG tiles unchanged, builds the lower pyramid levels, verifies the
result with OpenSlide, and writes a summary CSV.

## Locate the script

Use the first of these that exists:
1. `$DICOM2SVS` (environment variable holding the full path to `dicom2svs.py`)
2. `./dicom2svs.py` (when working inside the cloned repo)
3. `~/dicom2svs/dicom2svs.py`

If none exists, ask the user where they cloned the repo (https://github.com/GJiananChen/dicom2svs).
It needs `uv`. If `uv --version` fails, tell the user to install it
(`curl -LsSf https://astral.sh/uv/install.sh | sh` or `brew install uv`) instead of pip-installing packages.

## Workflow

1. **Get the inputs.** You need one or more source folders (scanned recursively) and an
   output folder. If the user gives no output folder, suggest `<first source folder>/../svs_output`
   and confirm it with them.
2. **Always dry-run first:**
   ```bash
   uv run "$SCRIPT" SRC1 [SRC2 ...] -o OUT --dry-run
   ```
   Show the user how many slides were found and which ones pass or fail the pre-checks.
   Explain each failure in plain words:
   - `Full resolution … DERIVED/RESAMPLED`: only lower-resolution levels were exported, so the full-resolution file is missing.
   - `Objective 40x … 20x` or `40x pixel size`: the scan is not 40x. The tool converts 40x only, by design.
   - `Tile encoding: unsupported`: JPEG 2000 or another codec. Suggest re-exporting as JPEG baseline.
   - `Whole-slide image`: not pathology WSI (e.g. CT/MR). This is expected in mixed folders.
3. **Convert** once the user confirms. Run in the background for big batches, since each GB of DICOM takes roughly 20–60 s:
   ```bash
   uv run "$SCRIPT" SRC1 [SRC2 ...] -o OUT
   ```
   Add `--overwrite` only if the user asks to replace existing `.svs` files.
4. **Report** from the final summary table: counts of converted, skipped, failed and exists;
   the output paths; and the location of the summary CSV (`OUT/dicom2svs_summary_*.csv`).
   The post-conversion checks prove the output is 40x, full size, has the same MPP, and has
   pixels identical to the source. Say so only when those checks printed `[ok]`. Exit code 1
   means at least one slide failed; quote its reason from the "Details" section.

## Don'ts

- Don't `pip install` the dependencies. `uv run` reads them from the script header.
- Don't loosen the 40x checks or edit the script to force a conversion unless the user explicitly asks.
- Don't print or upload slide label images or DICOM patient fields. The tool doesn't copy labels.
