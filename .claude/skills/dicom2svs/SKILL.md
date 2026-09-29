---
name: dicom2svs
description: Convert DICOM whole-slide images (pathology WSI, e.g. Leica/Aperio GT450 or Sectra exports) into pyramidal Aperio .svs files, one per folder from its largest DICOM file, flagging metadata that doesn't look like 40x full resolution for review. Use when the user wants to convert DICOM slides/folders to SVS, find the full-resolution DICOM in a slide folder, or batch-convert WSI DICOM for QuPath/OpenSlide.
---

# dicom2svs

Converts the largest DICOM file of each folder (normally the 40x full-resolution scan) into
one Aperio `.svs`. It copies the full-resolution JPEG tiles unchanged, builds the lower
pyramid levels, verifies the result with OpenSlide, and writes a summary CSV.

Metadata problems (no objective power, pixel size outside the 40x range, missing ImageType, …)
are **warnings**: the slide is still converted and the warning is recorded in the CSV
(`needs_review` = yes). Only unreadable files and unsupported tile encodings are skipped.
`--strict` skips warned slides instead.

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
   Show the user how many folders were found, which will convert cleanly, which will
   convert with warnings, and which will be skipped. Explain them in plain words:
   - `[WARN] Objective 40x … not recorded`: common in anonymised exports; the magnification is inferred from the pixel size (`40x*`).
   - `[WARN] Objective 40x … 20x` or `[WARN] 40x pixel size`: probably a real non-40x scan. It is still converted; offer `--strict` if the user wants 40x only.
   - `[WARN] Full resolution … DERIVED/RESAMPLED`: the largest file is a lower level, so the full-resolution file is probably missing from the export.
   - `[WARN] Magnification: .svs has no magnification or MPP`: neither was recorded; the user must set the pixel size in their viewer.
   - `[FAIL] Tile encoding: unsupported` (skipped): JPEG 2000 or another codec. Suggest re-exporting as JPEG baseline.
   - `[FAIL] Readable` (skipped): damaged or incomplete file.
3. **Convert** once the user confirms. Run in the background for big batches, since each GB of DICOM takes roughly 20–60 s:
   ```bash
   uv run "$SCRIPT" SRC1 [SRC2 ...] -o OUT
   ```
   Add `--overwrite` only if the user asks to replace existing `.svs` files, and `--strict`
   only if they want slides with warnings skipped.
4. **Report** from the final summary table: counts of converted, skipped, failed and exists;
   how many need review; the output paths; and the location of the summary CSV
   (`OUT/dicom2svs_summary_*.csv`, filter on `needs_review`). List the slides with warnings
   and their reasons. The post-conversion checks prove the output is full size, has the
   magnification and MPP that were written, and has pixels identical to the source. Say so only when those checks printed `[ok]`. Exit code 1
   means at least one slide failed; quote its reason from the "Details" section.

## Don'ts

- Don't `pip install` the dependencies. `uv run` reads them from the script header.
- Don't describe a slide with warnings as a verified 40x slide; point the user to its warnings.
- Don't edit the script to change which checks warn or skip unless the user explicitly asks.
- Don't print or upload slide label images or DICOM patient fields. The tool doesn't copy labels.
