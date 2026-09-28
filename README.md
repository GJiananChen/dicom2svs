# dicom2svs

Convert DICOM whole-slide images (WSI) into pyramidal Aperio `.svs` files, which open in
QuPath, OpenSlide, ASAP and most pathology AI pipelines.

Point it at one or more folders. For each slide it finds, it:

1. **Finds the largest DICOM file**, which is the full-resolution scan.
2. **Checks that it is 40x full resolution** before converting (details below).
3. **Copies the JPEG tiles into the `.svs` without re-compressing them**, so the 40x
   pixels are identical to the DICOM source.
4. **Builds the lower-resolution pyramid levels** (1/4, 1/16, 1/64 …) and a thumbnail.
5. **Checks the result** by opening it with OpenSlide and comparing it against the source.
6. **Prints a summary table** and saves it as a CSV in the output folder.

A progress bar shows each stage as it runs.

```
[1/1] ANON77746L14E
  largest file: /data/ANON77746L14E/ANON77746L14E_1_1.dcm (1.1 GB)
  Pre-conversion checks:
    [ok]   Whole-slide image: VL Whole Slide Microscopy
    [ok]   Full resolution: ImageType ORIGINAL/PRIMARY/VOLUME/NONE
    [ok]   Largest level: 156,123 x 76,339 px is the largest of 4 pyramid file(s)
    [ok]   Objective 40x: objective lens power 40x
    [ok]   40x pixel size: 0.2625 µm/px (40x range 0.2-0.3)
    [ok]   Tile encoding: 182,390 JPEG tiles of 256x256, single focal plane
  Post-conversion checks (OpenSlide):
    [ok]   Opens as SVS: OpenSlide vendor 'aperio'
    [ok]   Full resolution kept: 156,123 x 76,339 px
    [ok]   Reads as 40x: objective power 40
    [ok]   MPP kept: 0.2625 µm/px
    [ok]   Pixels identical: centre tile max difference 0 (tiles copied without re-compression)
    [ok]   Pyramid: 5 levels, downsamples 1, 3.99996, 15.9984, 63.987, 255.627
  -> /data/svs/ANON77746L14E.svs (1.2 GB, 20s)

==============================================================================
CONVERSION SUMMARY
==============================================================================
Slide          Status     Pixels          Mag  MPP     DICOM   SVS     Time
-------------  ---------  --------------  ---  ------  ------  ------  ----
ANON77746L14E  converted  156,123x76,339  40x  0.2625  1.1 GB  1.2 GB  20s
```

The example above is a Leica Aperio GT450 DX slide (1.1 GB) exported from Sectra. It
converted in about 20 seconds on an Apple Silicon Mac.

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

## Usage

```bash
# Check first: scan the folders and run the checks, but don't write anything
uv run dicom2svs.py /path/to/dicom_folder_1 /path/to/dicom_folder_2 -o /path/to/svs_out --dry-run

# Convert
uv run dicom2svs.py /path/to/dicom_folder_1 /path/to/dicom_folder_2 -o /path/to/svs_out
```

| Option | Default | Meaning |
|---|---|---|
| `DIR [DIR ...]` | (required) | Folders to scan. Subfolders are included. |
| `-o, --output` | (required) | Folder to write `.svs` files and the summary CSV into. It is created if missing. |
| `--dry-run` | off | Scan and check only. Nothing is written. |
| `--overwrite` | off | Replace `.svs` files that already exist. Without it, they are skipped. |
| `--quality` | 90 | JPEG quality for the **generated** lower levels. The 40x level is always copied unchanged. |
| `--workers` | all cores | Threads used to build the pyramid levels. |

The exit code is `0` when nothing failed and `1` when any slide failed conversion or
verification. Slides skipped by the pre-checks don't count as failures.

### How slides are found and named

* The tool treats every file with the DICOM `DICM` signature as DICOM, whatever its
  extension. Other files such as `.import`, `.zip` or `.jpg` are ignored.
* DICOM files in the **same folder** with the **same SeriesInstanceUID** count as one slide.
  The largest file by size is the one checked and converted.
* The output is named after the folder, e.g. `ANON77746L14E/` becomes `ANON77746L14E.svs`.
  If a folder holds several slides, the file name is appended to keep them apart.
* The file is written to `NAME.svs.partial` first and renamed once complete, so an
  interrupted run never leaves a half-written `.svs`.

## The 40x full-resolution check

The tool converts a slide only if its largest DICOM file passes **all** of these checks.
Any failure skips the slide and gives the reason in the summary.

| Check | Rule |
|---|---|
| Whole-slide image | SOP Class is *VL Whole Slide Microscopy Image Storage* |
| Full resolution | `ImageType` is `ORIGINAL\…\VOLUME\NONE`, not `DERIVED`/`RESAMPLED`, a thumbnail, a label or an overview |
| Largest level | No other file of the slide has a larger pixel matrix |
| Objective 40x | `OpticalPathSequence › ObjectiveLensPower` = 40 |
| 40x pixel size | Pixel spacing is 0.20–0.30 µm/px (40x scanners are typically 0.22–0.275) |
| Tile encoding | JPEG baseline, `TILED_FULL`, a single focal plane and optical path, 8-bit RGB |

After conversion, OpenSlide reopens the `.svs`. It must report vendor `aperio`, objective
power 40, the same pixel dimensions and MPP as the source, and a centre tile that decodes
**identically** to the matching DICOM frame.

## Output format

The output follows the standard Aperio SVS layout, the same page order scanners write:

| TIFF page | Content |
|---|---|
| 0 | Full-resolution (40x) level: the original JPEG tiles, plus the ICC colour profile and an Aperio description (`AppMag = 40`, `MPP = …`) |
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
  scanner model, MPP and magnification. No other patient fields are copied.

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
Claude does a dry run first, shows you the check results, and converts only after you
confirm.

## Troubleshooting

| Symptom | Fix |
|---|---|
| `uv: command not found` | Install uv (Setup, step 1), then open a new terminal. |
| `Tile encoding: unsupported: JPEG 2000 …` | This export uses a codec the tool can't copy. Re-export from the PACS as JPEG baseline if you can. |
| `Objective 40x: objective lens power 20x` | This is a 20x scan. The tool converts 40x scans only, by design. |
| `Full resolution: ImageType DERIVED/…/RESAMPLED` | The folder has only lower-resolution levels. The full-resolution file is missing from the export. |
| `… already exists (use --overwrite …)` | Delete the old `.svs` or rerun with `--overwrite`. |
| Slow on a network drive | Copy the DICOM folder to a local disk first. The tool reads the full-resolution file twice. |

## License

[MIT](LICENSE) © 2026 Jianan Chen
