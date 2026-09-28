#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.10"
# dependencies = [
#     "pydicom>=3.0",
#     "tifffile>=2024.1.30",
#     "imagecodecs",
#     "pillow>=10",
#     "numpy",
#     "tqdm",
#     "openslide-python",
#     "openslide-bin",
# ]
# ///
"""
dicom2svs — find the full-resolution DICOM whole-slide image in one or more
directory trees and convert it to an Aperio-style pyramidal .svs file.

For every slide found (DICOM files grouped by folder and SeriesInstanceUID) the tool:
  1. picks the largest DICOM file of the slide,
  2. checks that it is a 40x, full-resolution (ORIGINAL/VOLUME) image,
  3. copies its JPEG tiles into the .svs unchanged (no re-compression),
  4. builds the lower-resolution pyramid levels and a thumbnail from it,
  5. re-opens the result with OpenSlide to confirm it reads as a 40x slide
     with the same dimensions and pixels as the source.

Usage:
  uv run dicom2svs.py DIR [DIR ...] -o TARGET_DIR [--dry-run] [--overwrite]
"""

from __future__ import annotations

import argparse
import csv
import io
import math
import os
import struct
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path

import numpy as np
import pydicom
import tifffile
from PIL import Image
from tqdm import tqdm

WSI_SOP_CLASS = "1.2.840.10008.5.1.4.1.1.77.1.6"  # VL Whole Slide Microscopy Image Storage
JPEG_BASELINE = "1.2.840.10008.1.2.4.50"
EXPECTED_OBJECTIVE = 40.0
MPP_RANGE_40X = (0.20, 0.30)  # µm/px; 40x scanners are ~0.22-0.275 µm/px
PYRAMID_FACTOR = 4
PYRAMID_TILE = 256
PYRAMID_MIN_SIZE = 1024  # stop adding levels once the level fits in this many px
THUMBNAIL_SIZE = 1024
CLASSIC_TIFF_LIMIT = 3.8 * 1024**3  # switch to BigTIFF above this estimated size
WHITE = (255, 255, 255)


# --------------------------------------------------------------------------- #
# Scanning
# --------------------------------------------------------------------------- #

@dataclass
class DicomFile:
    path: Path
    size: int
    series_uid: str
    image_type: tuple[str, ...]
    width: int
    height: int


@dataclass
class Slide:
    series_uid: str
    files: list[DicomFile]
    name: str = ""
    largest: DicomFile | None = None


def is_dicom(path: Path) -> bool:
    try:
        with open(path, "rb") as f:
            f.seek(128)
            return f.read(4) == b"DICM"
    except OSError:
        return False


def scan(dirs: list[Path]) -> list[Slide]:
    """Find every DICOM file under `dirs` and group them into slides by folder and series."""
    candidates = []
    for d in dirs:
        for root, subdirs, names in os.walk(d):
            subdirs[:] = sorted(s for s in subdirs if not s.startswith("."))
            candidates += [Path(root) / n for n in sorted(names) if not n.startswith(".")]

    slides: dict[tuple[Path, str], Slide] = {}
    for path in tqdm(candidates, desc="Scanning", unit="file", leave=False):
        if not is_dicom(path):
            continue
        try:
            ds = pydicom.dcmread(
                path, stop_before_pixels=True,
                specific_tags=["SeriesInstanceUID", "ImageType",
                               "TotalPixelMatrixColumns", "TotalPixelMatrixRows",
                               "Columns", "Rows"],
            )
        except Exception:
            continue
        uid = str(ds.get("SeriesInstanceUID", "") or f"no-series:{path.parent}")
        f = DicomFile(
            path=path,
            size=path.stat().st_size,
            series_uid=uid,
            image_type=tuple(ds.get("ImageType", ())),
            width=int(ds.get("TotalPixelMatrixColumns", ds.get("Columns", 0)) or 0),
            height=int(ds.get("TotalPixelMatrixRows", ds.get("Rows", 0)) or 0),
        )
        slides.setdefault((path.parent, uid), Slide(series_uid=uid, files=[])).files.append(f)

    result = list(slides.values())
    for s in result:
        s.largest = max(s.files, key=lambda f: f.size)
    assign_names(result)
    return result


def assign_names(slides: list[Slide]) -> None:
    """Name each slide after its folder; disambiguate folders holding several slides."""
    by_dir: dict[Path, list[Slide]] = {}
    for s in slides:
        by_dir.setdefault(s.largest.path.parent, []).append(s)
    used: set[str] = set()
    for folder, group in by_dir.items():
        for s in group:
            name = folder.name
            if len(group) > 1 or name in used:
                name = f"{name}_{s.largest.path.stem}"
            while name in used:
                name += "_x"
            used.add(name)
            s.name = name


# --------------------------------------------------------------------------- #
# Pre-conversion checks
# --------------------------------------------------------------------------- #

@dataclass
class SourceInfo:
    path: Path
    width: int = 0
    height: int = 0
    tile_w: int = 0
    tile_h: int = 0
    frames: int = 0
    mpp: float = 0.0
    objective: float | None = None
    manufacturer: str = ""
    model: str = ""
    icc: bytes | None = None
    pixel_data_offset: int = 0
    jpeg_colorspace: str = ""
    subsampling: tuple[int, int] | None = None


@dataclass
class CheckResult:
    passed: list[str] = field(default_factory=list)
    failed: list[str] = field(default_factory=list)

    def check(self, ok: bool, label: str, detail: str) -> bool:
        (self.passed if ok else self.failed).append(f"{label}: {detail}")
        return ok

    @property
    def ok(self) -> bool:
        return not self.failed


def jpeg_colour_info(data: bytes) -> tuple[str, tuple[int, int] | None]:
    """Return ('ycbcr'|'rgb'|'unknown', chroma subsampling) of a baseline JPEG."""
    jfif = adobe_transform = None
    comps = []
    i = 2
    while i + 4 <= len(data) and data[i] == 0xFF:
        marker = data[i + 1]
        length = struct.unpack(">H", data[i + 2:i + 4])[0]
        seg = data[i + 4:i + 2 + length]
        if marker == 0xE0 and seg[:5] == b"JFIF\x00":
            jfif = True
        elif marker == 0xEE and seg[:5] == b"Adobe":
            adobe_transform = seg[11]
        elif marker in (0xC0, 0xC1):
            n = seg[5]
            comps = [(seg[6 + 3 * k], seg[7 + 3 * k] >> 4, seg[7 + 3 * k] & 15) for k in range(n)]
        elif marker == 0xDA:
            break
        i += 2 + length
    if len(comps) != 3:
        return "unknown", None
    if adobe_transform == 0 or [c[0] for c in comps] == [ord("R"), ord("G"), ord("B")]:
        return "rgb", None
    if jfif or adobe_transform == 1 or [c[0] for c in comps] == [1, 2, 3]:
        (_, h0, v0), (_, h1, v1), (_, h2, v2) = comps
        if h1 == h2 == v1 == v2 == 1:
            return "ycbcr", (h0, v0)
    return "unknown", None


def pixel_data_offset(path: Path) -> tuple[pydicom.Dataset, int]:
    with open(path, "rb") as f:
        ds = pydicom.dcmread(f, stop_before_pixels=True)
        offset = f.tell()
        f.seek(offset)
        tag = f.read(4)
    if tag != b"\xe0\x7f\x10\x00":
        raise ValueError("pixel data element not found where expected")
    return ds, offset


def first_frame(path: Path, offset: int) -> bytes:
    with open(path, "rb") as f:
        f.seek(offset + 12)  # skip (7FE0,0010) OB 0000 FFFFFFFF
        return next(pydicom.encaps.generate_frames(f, number_of_frames=1))


def check_source(slide: Slide) -> tuple[CheckResult, SourceInfo]:
    """Confirm the largest file of the slide is a 40x full-resolution WSI we can convert."""
    f = slide.largest
    res = CheckResult()
    info = SourceInfo(path=f.path)
    try:
        ds, info.pixel_data_offset = pixel_data_offset(f.path)
    except Exception as e:
        res.check(False, "Readable", f"cannot read DICOM header ({e})")
        return res, info

    res.check(ds.get("SOPClassUID") == WSI_SOP_CLASS, "Whole-slide image",
              "VL Whole Slide Microscopy" if ds.get("SOPClassUID") == WSI_SOP_CLASS
              else f"SOP class is {getattr(ds.get('SOPClassUID'), 'name', 'unknown')}")

    image_type = [str(t).upper() for t in ds.get("ImageType", [])]
    full_res = (len(image_type) >= 3 and image_type[0] == "ORIGINAL" and image_type[2] == "VOLUME"
                and (len(image_type) < 4 or image_type[3] == "NONE"))
    res.check(full_res, "Full resolution",
              f"ImageType {'/'.join(image_type) or 'missing'}"
              + ("" if full_res else " (expected ORIGINAL/.../VOLUME/NONE)"))

    info.width = int(ds.get("TotalPixelMatrixColumns", 0) or 0)
    info.height = int(ds.get("TotalPixelMatrixRows", 0) or 0)
    siblings = [s for s in slide.files if s is not f and "VOLUME" in s.image_type]
    bigger = [s for s in siblings if s.width * s.height > info.width * info.height]
    res.check(not bigger, "Largest level",
              f"{info.width:,} x {info.height:,} px is the largest of {len(siblings) + 1} "
              f"pyramid file(s)" if not bigger else f"{bigger[0].path.name} has more pixels")

    try:
        info.objective = float(ds.OpticalPathSequence[0].ObjectiveLensPower)
    except Exception:
        info.objective = None
    res.check(info.objective == EXPECTED_OBJECTIVE, "Objective 40x",
              f"objective lens power {info.objective:g}x" if info.objective is not None
              else "objective lens power not recorded")

    try:
        spacing = ds.SharedFunctionalGroupsSequence[0].PixelMeasuresSequence[0].PixelSpacing
        info.mpp = float(spacing[1]) * 1000  # mm -> µm; PixelSpacing is (row, column)
        mpp_y = float(spacing[0]) * 1000
    except Exception:
        info.mpp = mpp_y = 0.0
    lo, hi = MPP_RANGE_40X
    res.check(lo <= info.mpp <= hi and abs(info.mpp - mpp_y) < 1e-3, "40x pixel size",
              f"{info.mpp:.4f} µm/px (40x range {lo}-{hi})" if info.mpp
              else "pixel spacing not recorded")

    syntax = ds.file_meta.get("TransferSyntaxUID")
    info.tile_w, info.tile_h = int(ds.get("Columns", 0)), int(ds.get("Rows", 0))
    info.frames = int(ds.get("NumberOfFrames", 1))
    expected_frames = (math.ceil(info.width / max(info.tile_w, 1))
                       * math.ceil(info.height / max(info.tile_h, 1)))
    layout_ok = (
        syntax == JPEG_BASELINE
        and ds.get("DimensionOrganizationType", "TILED_FULL") == "TILED_FULL"
        and int(ds.get("TotalPixelMatrixFocalPlanes", 1)) == 1
        and int(ds.get("NumberOfOpticalPaths", 1)) == 1
        and int(ds.get("SamplesPerPixel", 0)) == 3
        and int(ds.get("BitsAllocated", 0)) == 8
        and info.frames == expected_frames
    )
    if layout_ok:
        info.jpeg_colorspace, info.subsampling = jpeg_colour_info(
            first_frame(f.path, info.pixel_data_offset))
        layout_ok = info.jpeg_colorspace != "unknown"
    res.check(layout_ok, "Tile encoding",
              f"{info.frames:,} JPEG tiles of {info.tile_w}x{info.tile_h}, single focal plane"
              if layout_ok else
              f"unsupported: {getattr(syntax, 'name', syntax)}, "
              f"{ds.get('DimensionOrganizationType', '?')}, {info.frames} frames "
              f"(this tool copies TILED_FULL JPEG-baseline tiles)")

    info.manufacturer = str(ds.get("Manufacturer", ""))
    info.model = str(ds.get("ManufacturerModelName", ""))
    try:
        info.icc = bytes(ds.OpticalPathSequence[0].ICCProfile)
    except Exception:
        info.icc = None
    return res, info


# --------------------------------------------------------------------------- #
# Conversion
# --------------------------------------------------------------------------- #

class TileSource:
    """Random access to the encapsulated JPEG frames of a DICOM file."""

    def __init__(self, info: SourceInfo):
        self.info = info
        self.fd = os.open(info.path, os.O_RDONLY)
        self.index = self._build_index()
        self.cols = math.ceil(info.width / info.tile_w)
        self.rows = math.ceil(info.height / info.tile_h)

    def _build_index(self) -> list[tuple[int, int]]:
        pos = self.info.pixel_data_offset + 12
        index = []
        end = os.fstat(self.fd).st_size
        first = True
        while pos + 8 <= end:
            group, elem, length = struct.unpack("<HHI", os.pread(self.fd, 8, pos))
            pos += 8
            if (group, elem) == (0xFFFE, 0xE0DD):  # sequence delimiter
                break
            if (group, elem) != (0xFFFE, 0xE000):
                raise ValueError("corrupt encapsulated pixel data")
            if not first:  # the first item is the basic offset table
                index.append((pos, length))
            first = False
            pos += length
        if len(index) != self.info.frames:
            raise ValueError(f"found {len(index)} fragments for {self.info.frames} frames "
                             "(multi-fragment frames are not supported)")
        return index

    def tile(self, row: int, col: int) -> bytes:
        offset, length = self.index[row * self.cols + col]
        return os.pread(self.fd, length, offset)

    def close(self) -> None:
        os.close(self.fd)


@dataclass
class Level:
    width: int
    height: int
    tile_w: int
    tile_h: int
    tiles: list[bytes]  # row-major encoded JPEG tiles

    @property
    def cols(self) -> int:
        return math.ceil(self.width / self.tile_w)

    @property
    def rows(self) -> int:
        return math.ceil(self.height / self.tile_h)

    def tile(self, row: int, col: int) -> bytes:
        return self.tiles[row * self.cols + col]


def decode(data: bytes) -> Image.Image:
    return Image.open(io.BytesIO(data)).convert("RGB")


def read_region(src, width, height, tile_w, tile_h, x0, y0, w, h) -> Image.Image:
    """Assemble a w x h region at (x0, y0) from a tiled source; outside the image is white."""
    canvas = Image.new("RGB", (w, h), WHITE)
    x1, y1 = min(x0 + w, width), min(y0 + h, height)
    for r in range(y0 // tile_h, (y1 - 1) // tile_h + 1):
        for c in range(x0 // tile_w, (x1 - 1) // tile_w + 1):
            canvas.paste(decode(src.tile(r, c)), (c * tile_w - x0, r * tile_h - y0))
    if x1 - x0 < w:  # blank out the scanner's edge-tile padding
        canvas.paste(WHITE, (x1 - x0, 0, w, h))
    if y1 - y0 < h:
        canvas.paste(WHITE, (0, y1 - y0, w, h))
    return canvas


def build_level(src, width, height, tile_w, tile_h, quality, subsampling, workers, bar_desc) -> Level:
    """Downsample a tiled source by PYRAMID_FACTOR into a new level of JPEG tiles."""
    f, t = PYRAMID_FACTOR, PYRAMID_TILE
    level = Level(math.ceil(width / f), math.ceil(height / f), t, t, [])
    jpeg_sub = {(1, 1): 0, (2, 1): 1, (2, 2): 2}.get(subsampling, 2)

    def make(i: int) -> bytes:
        r, c = divmod(i, level.cols)
        region = read_region(src, width, height, tile_w, tile_h, c * t * f, r * t * f, t * f, t * f)
        out = io.BytesIO()
        region.reduce(f).save(out, "JPEG", quality=quality, subsampling=jpeg_sub)
        return out.getvalue()

    n = level.cols * level.rows
    with ThreadPoolExecutor(workers) as pool, \
            tqdm(total=n, desc=bar_desc, unit="tile", leave=False, dynamic_ncols=True) as bar:
        for data in pool.map(make, range(n), chunksize=16):
            level.tiles.append(data)
            bar.update()
    return level


def aperio_description(info: SourceInfo, name: str, quality: int, level: Level | None) -> str:
    w, h = info.width, info.height
    head = f"Aperio Image Library v12.0.15 \r\n{w}x{h} [0,0 {w}x{h}] ({info.tile_w}x{info.tile_h})"
    if level is not None:
        return f"{head} -> {level.width}x{level.height} JPEG/RGB Q={quality}"
    fields = [
        f"{head} JPEG/RGB Q={quality}",
        f"AppMag = {info.objective:g}",
        f"MPP = {info.mpp:.6f}",
        f"Filename = {name}",
        f"Source = DICOM {info.path.name}",
        f"Scanner = {info.manufacturer} {info.model}".strip(),
        f"Converted = dicom2svs {datetime.now():%Y-%m-%d %H:%M}",
    ]
    return "|".join(fields)


def convert(slide: Slide, info: SourceInfo, out_path: Path, quality: int, workers: int,
            label: str) -> None:
    src = TileSource(info)
    try:
        # Reduced levels are built first so the thumbnail can be placed right after
        # the base image, matching the page order of scanner-written .svs files.
        levels: list[Level] = []
        prev = (src, info.width, info.height, info.tile_w, info.tile_h)
        while max(prev[1], prev[2]) > PYRAMID_MIN_SIZE:
            n = len(levels) + 1
            lvl = build_level(*prev, quality, info.subsampling, workers,
                              f"{label} level {n} (1/{PYRAMID_FACTOR ** n})")
            levels.append(lvl)
            prev = (lvl, lvl.width, lvl.height, lvl.tile_w, lvl.tile_h)

        smallest = levels[-1] if levels else None
        if smallest:
            thumb = read_region(smallest, smallest.width, smallest.height,
                                smallest.tile_w, smallest.tile_h, 0, 0,
                                smallest.width, smallest.height)
        else:
            thumb = read_region(src, info.width, info.height, info.tile_w, info.tile_h,
                                0, 0, info.width, info.height)
        thumb.thumbnail((THUMBNAIL_SIZE, THUMBNAIL_SIZE))

        estimate = sum(l for _, l in src.index) + sum(len(t) for l in levels for t in l.tiles)
        colour = dict(photometric="ycbcr", subsampling=info.subsampling) \
            if info.jpeg_colorspace == "ycbcr" else dict(photometric="rgb")
        extratags = [(34675, 7, len(info.icc), info.icc, True)] if info.icc else []

        partial = out_path.with_name(out_path.name + ".partial")
        with tifffile.TiffWriter(partial, bigtiff=estimate > CLASSIC_TIFF_LIMIT) as tw, \
                tqdm(total=len(src.index), desc=f"{label} level 0 (full res)", unit="tile",
                     leave=False, dynamic_ncols=True) as bar:

            def base_tiles():
                for r in range(src.rows):
                    for c in range(src.cols):
                        yield src.tile(r, c)
                        bar.update()

            tw.write(base_tiles(), shape=(info.height, info.width, 3), dtype="uint8",
                     tile=(info.tile_h, info.tile_w), compression="jpeg", **colour,
                     description=aperio_description(info, out_path.name, quality, None),
                     extratags=extratags, metadata=None)
            tw.write(np.asarray(thumb), compression="zlib", photometric="rgb",
                     rowsperstrip=16, metadata=None,
                     description=f"Aperio Image Library v12.0.15 \r\n{info.width}x{info.height}"
                                 f" -> {thumb.width}x{thumb.height} - thumbnail")
            for lvl in levels:
                tw.write(iter(lvl.tiles), shape=(lvl.height, lvl.width, 3), dtype="uint8",
                         tile=(lvl.tile_h, lvl.tile_w), compression="jpeg", **colour,
                         description=aperio_description(info, out_path.name, quality, lvl),
                         metadata=None)
        os.replace(partial, out_path)
    except BaseException:
        partial = out_path.with_name(out_path.name + ".partial")
        if partial.exists():
            partial.unlink()
        raise
    finally:
        src.close()


# --------------------------------------------------------------------------- #
# Post-conversion verification
# --------------------------------------------------------------------------- #

def verify(out_path: Path, info: SourceInfo) -> CheckResult:
    """Re-open the .svs with OpenSlide and compare it against the DICOM source."""
    import openslide

    res = CheckResult()
    slide = openslide.OpenSlide(str(out_path))
    try:
        p = slide.properties
        res.check(p.get("openslide.vendor") == "aperio", "Opens as SVS",
                  f"OpenSlide vendor '{p.get('openslide.vendor')}'")
        res.check(slide.dimensions == (info.width, info.height), "Full resolution kept",
                  f"{slide.dimensions[0]:,} x {slide.dimensions[1]:,} px")
        res.check(p.get("openslide.objective-power") == f"{info.objective:g}" == "40",
                  "Reads as 40x", f"objective power {p.get('openslide.objective-power')}")
        mpp = float(p.get("openslide.mpp-x", 0))
        res.check(abs(mpp - info.mpp) < 1e-4, "MPP kept", f"{mpp:.4f} µm/px")

        # Pixel check: the centre tile must decode identically to the DICOM frame.
        src = TileSource(info)
        try:
            r, c = src.rows // 2, src.cols // 2
            expected = np.asarray(decode(src.tile(r, c)), dtype=np.int16)
        finally:
            src.close()
        got = np.asarray(slide.read_region((c * info.tile_w, r * info.tile_h), 0,
                                           (info.tile_w, info.tile_h)).convert("RGB"),
                         dtype=np.int16)
        diff = int(np.abs(expected - got).max())
        res.check(diff <= 2, "Pixels identical",
                  f"centre tile max difference {diff} (tiles copied without re-compression)")
        res.check(slide.level_count > 1 or max(slide.dimensions) <= PYRAMID_MIN_SIZE,
                  "Pyramid", f"{slide.level_count} levels, downsamples "
                  + ", ".join(f"{d:g}" for d in slide.level_downsamples))
    finally:
        slide.close()
    return res


# --------------------------------------------------------------------------- #
# Reporting
# --------------------------------------------------------------------------- #

@dataclass
class Outcome:
    slide: Slide
    status: str  # converted | skipped | failed | would convert | exists
    info: SourceInfo | None = None
    output: Path | None = None
    seconds: float = 0.0
    notes: list[str] = field(default_factory=list)


def human(n: float) -> str:
    for unit in ("B", "KB", "MB", "GB"):
        if n < 1024 or unit == "GB":
            return f"{n:.1f} {unit}" if unit != "B" else f"{n:.0f} B"
        n /= 1024
    return f"{n:.1f} TB"


def print_checks(title: str, res: CheckResult) -> None:
    tqdm.write(f"  {title}")
    for line in res.passed:
        tqdm.write(f"    [ok]   {line}")
    for line in res.failed:
        tqdm.write(f"    [FAIL] {line}")


def print_summary(outcomes: list[Outcome], target: Path, elapsed: float) -> Path | None:
    print("\n" + "=" * 78)
    print("CONVERSION SUMMARY")
    print("=" * 78)
    if not outcomes:
        print("No DICOM slides found.")
        return None
    rows = []
    for o in outcomes:
        i = o.info
        rows.append([
            o.slide.name,
            o.status,
            f"{i.width:,}x{i.height:,}" if i and i.width else "-",
            f"{i.objective:g}x" if i and i.objective else "-",
            f"{i.mpp:.4f}" if i and i.mpp else "-",
            human(o.slide.largest.size),
            human(o.output.stat().st_size) if o.output and o.output.exists() else "-",
            f"{o.seconds:.0f}s" if o.seconds else "-",
        ])
    header = ["Slide", "Status", "Pixels", "Mag", "MPP", "DICOM", "SVS", "Time"]
    widths = [max(len(str(r[k])) for r in rows + [header]) for k in range(len(header))]
    fmt = "  ".join(f"{{:<{w}}}" for w in widths)
    print(fmt.format(*header))
    print(fmt.format(*("-" * w for w in widths)))
    for r in rows:
        print(fmt.format(*r))

    problems = [o for o in outcomes if o.notes]
    if problems:
        print("\nDetails:")
        for o in problems:
            for n in o.notes:
                print(f"  {o.slide.name}: {n}")

    counts: dict[str, int] = {}
    for o in outcomes:
        counts[o.status] = counts.get(o.status, 0) + 1
    print("\n" + ", ".join(f"{v} {k}" for k, v in counts.items())
          + f" — {len(outcomes)} slide(s) in {elapsed:.0f}s")

    if not target.exists():
        return None
    report = target / f"dicom2svs_summary_{datetime.now():%Y%m%d_%H%M%S}.csv"
    with open(report, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["slide", "status", "source_file", "output_file", "width", "height",
                    "objective", "mpp", "dicom_bytes", "svs_bytes", "seconds", "notes"])
        for o in outcomes:
            i = o.info
            w.writerow([o.slide.name, o.status, o.slide.largest.path,
                        o.output if o.output and o.output.exists() else "",
                        i.width if i else "", i.height if i else "",
                        i.objective if i and i.objective else "", f"{i.mpp:.6f}" if i else "",
                        o.slide.largest.size,
                        o.output.stat().st_size if o.output and o.output.exists() else "",
                        f"{o.seconds:.1f}", " | ".join(o.notes)])
    print(f"Summary saved to {report}")
    return report


# --------------------------------------------------------------------------- #
# Main
# --------------------------------------------------------------------------- #

def main() -> int:
    ap = argparse.ArgumentParser(
        description="Convert the full-resolution 40x DICOM whole-slide image in each "
                    "slide folder to a pyramidal Aperio .svs file.")
    ap.add_argument("dirs", nargs="+", type=Path, help="directories to scan (recursively)")
    ap.add_argument("-o", "--output", type=Path, required=True, help="target directory for .svs files")
    ap.add_argument("--dry-run", action="store_true", help="scan and check only, do not convert")
    ap.add_argument("--overwrite", action="store_true", help="replace existing .svs files")
    ap.add_argument("--quality", type=int, default=90,
                    help="JPEG quality for the generated lower-resolution levels (default 90); "
                         "full-resolution tiles are always copied unchanged")
    ap.add_argument("--workers", type=int, default=os.cpu_count() or 4,
                    help="threads for building pyramid levels (default: all cores)")
    args = ap.parse_args()

    for d in args.dirs:
        if not d.is_dir():
            ap.error(f"not a directory: {d}")
    target = args.output.expanduser().resolve()

    start = time.time()
    slides = scan([d.expanduser().resolve() for d in args.dirs])
    print(f"Found {len(slides)} slide(s) in {sum(len(s.files) for s in slides)} DICOM file(s)")
    if not args.dry_run:
        target.mkdir(parents=True, exist_ok=True)

    outcomes: list[Outcome] = []
    for k, slide in enumerate(slides, 1):
        label = f"[{k}/{len(slides)}] {slide.name}"
        f = slide.largest
        tqdm.write(f"\n{label}\n  largest file: {f.path} ({human(f.size)})")
        checks, info = check_source(slide)
        print_checks("Pre-conversion checks:", checks)
        out_path = target / f"{slide.name}.svs"
        o = Outcome(slide, "", info, out_path)
        outcomes.append(o)

        if not checks.ok:
            o.status, o.output = "skipped", None
            o.notes += [f"check failed - {x}" for x in checks.failed]
            continue
        if out_path.exists() and not args.overwrite:
            o.status = "exists"
            o.notes.append(f"{out_path} already exists (use --overwrite to replace)")
            continue
        if args.dry_run:
            o.status, o.output = "would convert", None
            continue

        t0 = time.time()
        try:
            convert(slide, info, out_path, args.quality, args.workers, label)
            post = verify(out_path, info)
        except KeyboardInterrupt:
            o.status, o.output = "failed", None
            o.notes.append("interrupted by user")
            break
        except Exception as e:
            o.status, o.output = "failed", None
            o.notes.append(f"{type(e).__name__}: {e}")
            tqdm.write(f"  FAILED: {e}")
            continue
        finally:
            o.seconds = time.time() - t0
        print_checks("Post-conversion checks (OpenSlide):", post)
        if post.ok:
            o.status = "converted"
        else:
            o.status = "failed"
            o.notes += [f"verification failed - {x}" for x in post.failed]
        tqdm.write(f"  -> {out_path} ({human(out_path.stat().st_size)}, {o.seconds:.0f}s)")

    print_summary(outcomes, target, time.time() - start)
    return 1 if any(o.status == "failed" for o in outcomes) else 0


if __name__ == "__main__":
    sys.exit(main())
