# dicom2svs

Convert DICOM whole-slide images (WSI) into pyramidal Aperio `.svs` files, which open in
QuPath, OpenSlide, ASAP and most pathology AI pipelines.

Point it at one or more folders. For **each folder** that contains DICOM files, it:

1. **Takes the largest DICOM file**, which is normally the full-resolution scan. Each
   folder gives exactly one `.svs`.
2. **Checks the metadata** against a 40x full-resolution scan (details below). Missing or
   unexpected metadata is **recorded as a warning** for review; the slide is still
   converted. Only files that can't be converted at all are skipped.
3. **Copies the JPEG tiles into the `.svs` without re-compressing them**, so the
   full-resolution pixels are identical to the DICOM source.
4. **Builds the lower-resolution pyramid levels** (1/4, 1/16, 1/64 …) and a thumbnail.
5. **Checks the result** by opening it with OpenSlide and comparing it against the source.
6. **Prints a summary table** and saves a summary spreadsheet (CSV) in the output folder,
   with a `needs_review` column and every warning listed.

A progress bar shows each stage as it runs.

```
[1/2] ANON77746L14E
  largest file: /data/ANON77746L14E/ANON77746L14E_1_1.dcm (1.1 GB)
  Pre-conversion checks:
    [ok]   Whole-slide image: VL Whole Slide Microscopy
    [ok]   Full resolution: ImageType ORIGINAL/PRIMARY/VOLUME/NONE
    [ok]   Largest level: 156,123 x 76,339 px is the largest of 6 DICOM file(s) in the folder
    [ok]   Objective 40x: objective lens power 40x
    [ok]   40x pixel size: 0.2625 µm/px (40x range 0.2-0.3)
    [ok]   Image size: 156,123 x 76,339 px
    [ok]   Tile encoding: 182,390 JPEG tiles of 256x256, single focal plane
  Post-conversion checks (OpenSlide):
    [ok]   Opens as SVS: OpenSlide vendor 'aperio'
    [ok]   Full resolution kept: 156,123 x 76,339 px
    [ok]   Magnification kept: objective power 40 (metadata)
    [ok]   MPP kept: 0.2625 µm/px
    [ok]   Pixels identical: centre tile max difference 0 (tiles copied without re-compression)
    [ok]   Pyramid: 5 levels, downsamples 1, 3.99996, 15.9984, 63.987, 255.627
  -> /data/svs/ANON77746L14E.svs (1.2 GB, 19s)

[2/2] CASE_0042
  largest file: /data/CASE_0042/slide.dcm (1.3 GB)
  Pre-conversion checks:
    ...
    [WARN] Objective 40x: objective lens power not recorded
    [WARN] Magnification: .svs labelled 40x, inferred from 0.2630 µm/px
  ...

==============================================================================
CONVERSION SUMMARY
==============================================================================
Slide          Status     Pixels          Mag   MPP     DICOM   SVS     Time  Review
-------------  ---------  --------------  ----  ------  ------  ------  ----  ------------
ANON77746L14E  converted  156,123x76,339  40x   0.2625  1.1 GB  1.2 GB  19s   -
CASE_0042      converted  160,512x80,120  40x*  0.2630  1.3 GB  1.4 GB  22s   2 warning(s)
* magnification inferred from pixel size (objective not recorded)

2 converted — 2 slide(s) in 41s; 1 need review (see the needs_review column)
```

The first slide above is a Leica Aperio GT450 DX slide (1.1 GB) exported from Sectra. It
converted in about 20 seconds on an Apple Silicon Mac. The second (`CASE_0042`) is an
illustrative example of how a slide
with missing metadata is reported.

---

## Setup (about 1 minute)

The only thing you need is [**uv**](https://docs.astral.sh/uv/). It is a single binary, and
the script declares its own dependencies, so there is no virtualenv, `pip install` or
conda environment to set up.

```bash
# 1. Install uv (skip this if `uv --version` already works)
curl -LsSf https://astral.sh/uv/install.sh | sh      # macOS / Linux
# or: brew install uv
# Windows: powershell -c "irm https://astral.sh/uv/install.ps1 | iex"

# 2. Get the tool
git clone https://github.com/GJiananChen/dicom2svs.git
cd dicom2svs

# 3. Check that it runs (the first run downloads dependencies, which takes ~10 s)
uv run dicom2svs.py --help
```

Dependencies come from the header of `dicom2svs.py`: pydicom, tifffile, imagecodecs,
Pillow, numpy, tqdm and OpenSlide. `openslide-bin` ships the OpenSlide library itself, so
you don't need Homebrew or apt packages. Python 3.10 or newer is needed; uv downloads one
if necessary.

**Works on Windows, Linux and macOS.** To check it on your machine without real data,
run `uv run tests/run_tests.py` (see [Development](#development)).

## Usage

```bash
# Check first: scan the folders and run the checks, but don't write anything
uv run dicom2svs.py /path/to/dicom_folder_1 /path/to/dicom_folder_2 -o /path/to/svs_out --dry-run

# Convert
uv run dicom2svs.py /path/to/dicom_folder_1 /path/to/dicom_folder_2 -o /path/to/svs_out
```

On Windows (PowerShell or Command Prompt), quote paths that contain spaces:

```powershell
uv run dicom2svs.py "D:\Slides\batch 1" "\\server\share\batch2" -o "D:\Slides\svs_out"
```

| Option | Default | Meaning |
|---|---|---|
| `DIR [DIR ...]` | (required) | Folders to scan. Subfolders are included. |
| `-o, --output` | (required) | Folder to write `.svs` files and the summary CSV into. It is created if missing. |
| `--dry-run` | off | Scan and check only. Nothing is written. |
| `--overwrite` | off | Replace `.svs` files that already exist. Without it, they are skipped. |
| `--strict` | off | Skip slides with metadata warnings instead of converting them. |
| `--quality` | 90 | JPEG quality for the **generated** lower levels. The full-resolution level is always copied unchanged. |
| `--workers` | all cores | Threads used to build the pyramid levels. |

The exit code is `0` when nothing failed and `1` when any slide failed conversion or
verification. Skipped slides and warnings don't count as failures.

### How slides are found and named

* The tool treats every file with the DICOM `DICM` signature as DICOM, whatever its
  extension. Other files such as `.import`, `.zip` or `.jpg` are ignored.
* **One folder = one slide = one `.svs`.** Only the largest DICOM file in each folder, by
  file size, is checked and converted. The other files (lower levels, label, overview)
  are left alone.
* The output is named after the folder, e.g. `ANON77746L14E/` becomes `ANON77746L14E.svs`.
  If two folders in different places have the same name, the parent folder's name is
  added to keep them apart.
* The file is written to `NAME.svs.partial` first and renamed once complete, so an
  interrupted run never leaves a half-written `.svs`.

## Checks: warnings and failures

Many exports have metadata that is missing or badly entered. So the 40x checks are
**warnings**: the slide is still converted, and each problem is listed in the summary
spreadsheet for you to review. Use `--strict` to skip warned slides instead.

| Check | Expected | If not met |
|---|---|---|
| Whole-slide image | SOP Class is *VL Whole Slide Microscopy Image Storage* | warning |
| Full resolution | `ImageType` is `ORIGINAL\…\VOLUME\NONE`, not `DERIVED`/`RESAMPLED`, a thumbnail, a label or an overview | warning |
| Largest level | No other file in the folder has a larger pixel matrix | warning |
| Objective 40x | `OpticalPathSequence › ObjectiveLensPower` = 40 | warning |
| 40x pixel size | Pixel spacing is 0.20–0.30 µm/px (40x scanners are typically 0.22–0.275) | warning |
| Magnification consistent | The objective matches the pixel size (≈10 ÷ µm/px) | warning |
| Readable | The DICOM header and pixel data can be read | **skipped** |
| Image size | `TotalPixelMatrixColumns/Rows` present | **skipped** |
| Tile encoding | JPEG baseline, `TILED_FULL`, a single focal plane and optical path, 8-bit RGB | **skipped** |

**Magnification written to the `.svs`:**
* If the objective power is recorded, that value is used.
* If it is missing but the pixel size is known, the magnification is inferred from the
  pixel size, e.g. 0.2625 µm/px gives 40x. It is marked `40x*` in the table and
  `inferred from MPP` in the spreadsheet.
* If both are missing, no magnification or MPP is written, and a warning says so.
  Set them in the viewer before measuring anything.

After conversion, OpenSlide reopens the `.svs`. It must report vendor `aperio`, the same
pixel dimensions, magnification and MPP that were written, and a centre tile that decodes
**identically** to the matching DICOM frame. If any of these fail, the slide is marked
`failed`.

### Summary spreadsheet

`dicom2svs_summary_<date>_<time>.csv` in the output folder has one row per folder:

| Column | Meaning |
|---|---|
| `slide`, `status` | Output name; `converted`, `skipped`, `failed`, `exists` or `would convert` |
| `needs_review` | `yes` if there are warnings or the slide was skipped or failed. Filter on this. |
| `warning_count`, `warnings` | Metadata problems, separated by ` \| ` |
| `errors` | Why a slide was skipped or failed |
| `source_folder`, `source_file`, `dicom_files_in_folder` | Where the slide came from |
| `output_file` | The `.svs` that was written |
| `width`, `height`, `mpp` | Full-resolution size and pixel size (µm/px) |
| `objective_recorded`, `svs_magnification`, `magnification_source` | Magnification from the DICOM, the one written to the `.svs`, and where it came from |
| `dicom_bytes`, `svs_bytes`, `seconds` | File sizes and conversion time |

**Note:** because warned slides are converted, a genuine 20x scan is converted too. It
shows up as `needs_review = yes` with the reason. Check that column after each batch, or
use `--strict`.

## Output format

The output follows the standard Aperio SVS layout, the same page order scanners write:

| TIFF page | Content |
|---|---|
| 0 | Full-resolution level: the original JPEG tiles, plus the ICC colour profile and an Aperio description (`AppMag = …`, `MPP = …` when known) |
| 1 | Thumbnail (at most 1024 px) |
| 2… | Pyramid levels, each 4x smaller than the one before, until the level fits in 1024 px |

It switches to BigTIFF automatically when the output would be larger than about 3.8 GB.

## Limitations

* **Only JPEG-baseline tiles are supported.** That covers Leica/Aperio GT450 exports and
  most Sectra exports. JPEG 2000, JPEG-LS and uncompressed DICOM WSI are skipped with a
  clear reason; they would need decoding and re-encoding.
* **One focal plane, brightfield only.** Z-stacks and multi-channel fluorescence are skipped.
* **Label and macro images are not copied**, since they often show identifiers.
  Only the scan itself goes into the `.svs`.
* **Some patient details are copied.** The `.svs` description stores the DICOM file name,
  the output (folder) name, the scanner model, MPP and magnification. No other patient
  fields are copied.

---

## Using it with Claude Code

This repo includes a Claude Code **skill** (`.claude/skills/dicom2svs/SKILL.md`) and a
`CLAUDE.md`. Claude therefore already knows how to run the tool, read the checks and
explain the summary.

**Option A: work inside the repo.** Run `claude` in the cloned folder and ask in plain
language:

> Convert the DICOM slides in ~/data/batch1 and ~/data/batch2 to svs in ~/data/svs

**Option B: use it from any project.** Install the skill for your user so it works
everywhere:

```bash
mkdir -p ~/.claude/skills
cp -r .claude/skills/dicom2svs ~/.claude/skills/
# The skill expects the script at ~/dicom2svs/dicom2svs.py by default. If you cloned
# somewhere else, set DICOM2SVS to the full path of dicom2svs.py, e.g. in ~/.zshrc:
export DICOM2SVS="$HOME/path/to/dicom2svs/dicom2svs.py"
```

Then type `/dicom2svs` in Claude Code, or just ask it to convert DICOM slides to SVS.
Claude does a dry run first, shows you the check results and warnings, and converts only
after you confirm.

## Troubleshooting

| Message | Meaning / fix |
|---|---|
| `uv: command not found` | Install uv (Setup, step 1), then open a new terminal. |
| `Tile encoding: unsupported: JPEG 2000 …` (skipped) | This export uses a codec the tool can't copy. Re-export from the PACS as JPEG baseline if you can. |
| `Readable: cannot read DICOM header` (skipped) | The file is damaged or not a complete DICOM file. Re-export it. |
| `[WARN] Objective 40x: objective lens power not recorded` | Common in anonymised exports. The magnification is inferred from the pixel size; check it in the spreadsheet. |
| `[WARN] Objective 40x: objective lens power 20x` | Probably a real 20x scan. It is converted as 20x; use `--strict` to skip these. |
| `[WARN] Full resolution: ImageType DERIVED/…/RESAMPLED` | The largest file in the folder is a lower-resolution level, so the full-resolution file is probably missing from the export. |
| `[WARN] Magnification: .svs has no magnification or MPP` | Neither was recorded. Set the pixel size in your viewer before measuring. |
| `… already exists (use --overwrite …)` | Delete the old `.svs` or rerun with `--overwrite`. |
| Slow on a network drive | Copy the DICOM folder to a local disk first. The tool reads the full-resolution file twice. |

## Development

`tests/run_tests.py` builds three small synthetic slides (no patient data), converts them
and checks the results: a clean 40x slide, one with missing metadata, and one with an
unsupported codec. Run it before pushing changes:

```bash
uv run tests/run_tests.py
```

`tests/make_test_slides.py OUT_DIR` writes the same synthetic slides if you want to try
the tool without real data.

## License

[MIT](LICENSE) © 2026 Jianan Chen
