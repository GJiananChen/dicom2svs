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

For every folder containing DICOM files the tool:
  1. picks the largest DICOM file in the folder (one .svs per folder),
  2. checks that it is a 40x, full-resolution (ORIGINAL/VOLUME) image; missing or
     unexpected metadata is recorded as a warning for review, not a reason to skip,
  3. copies its JPEG tiles into the .svs unchanged (no re-compression),
  4. builds the lower-resolution pyramid levels and a thumbnail from it,
  5. re-opens the result with OpenSlide to confirm it has the same dimensions,
     magnification, MPP and pixels as the source.

Usage:
  uv run dicom2svs.py DIR [DIR ...] -o TARGET_DIR [--dry-run] [--overwrite] [--strict]
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
import warnings
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path

import numpy as np
import pydicom
import tifffile
from PIL import Image
from tqdm import tqdm

# Badly entered values are reported by our own checks, so pydicom's warnings are noise.
pydicom.config.settings.reading_validation_mode = pydicom.config.IGNORE
warnings.filterwarnings("ignore", module="pydicom")

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
    image_type: tuple[str, ...]
    width: int
    height: int


@dataclass
class Slide:
    """All DICOM files of one folder; only the largest one is converted."""
    folder: Path
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
    """Find every DICOM file under `dirs`; each folder containing DICOM files is one slide."""
    candidates = []
    for d in dirs:
        for root, subdirs, names in os.walk(d):
            subdirs[:] = sorted(s for s in subdirs if not s.startswith("."))
            candidates += [Path(root) / n for n in sorted(names) if not n.startswith(".")]

    candidates = list(dict.fromkeys(candidates))  # a folder may sit under two given dirs

    slides: dict[Path, Slide] = {}
    for path in tqdm(candidates, desc="Scanning", unit="file", leave=False):
        if not is_dicom(path):
            continue
        try:
            ds = pydicom.dcmread(
                path, stop_before_pixels=True,
                specific_tags=["ImageType", "TotalPixelMatrixColumns", "TotalPixelMatrixRows",
                               "Columns", "Rows"],
            )
            image_type = tuple(str(t) for t in ds.get("ImageType", ()) or ())
            width = int(ds.get("TotalPixelMatrixColumns", ds.get("Columns", 0)) or 0)
            height = int(ds.get("TotalPixelMatrixRows", ds.get("Rows", 0)) or 0)
        except Exception:
            # Badly formed headers still count; the pre-checks report what is wrong.
            image_type, width, height = (), 0, 0
        f = DicomFile(path, path.stat().st_size, image_type, width, height)
        slides.setdefault(path.parent, Slide(folder=path.parent, files=[])).files.append(f)

    result = list(slides.values())
    for s in result:
        s.largest = max(s.files, key=lambda f: f.size)
    assign_names(result)
    return result


def assign_names(slides: list[Slide]) -> None:
    """Name each slide after its folder, adding the parent folder when names collide."""
    counts: dict[str, int] = {}
    for s in slides:
        counts[s.folder.name] = counts.get(s.folder.name, 0) + 1
    used: set[str] = set()
    for s in slides:
        name = s.folder.name
        if counts[name] > 1:
            name = f"{s.folder.parent.name}_{name}"
        base, n = name, 2
        while name in used:
            name, n = f"{base}_{n}", n + 1
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
    app_mag: float | None = None  # magnification written to the .svs
    mag_source: str = "unknown"   # metadata | inferred from MPP | unknown
    manufacturer: str = ""
    model: str = ""
    icc: bytes | None = None
    pixel_data_offset: int = 0
    jpeg_colorspace: str = ""
    subsampling: tuple[int, int] | None = None


@dataclass
class CheckResult:
    """`failed` stops the conversion; `warnings` are recorded for review only."""
    passed: list[str] = field(default_factory=list)
    failed: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)

    def check(self, ok: bool, label: str, detail: str, fatal: bool = True) -> bool:
        target = self.passed if ok else self.failed if fatal else self.warnings
        target.append(f"{label}: {detail}")
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


def safe(get, default=None):
    """Read one metadata value; badly entered or missing values give `default`."""
    try:
        value = get()
        return default if value is None or value == "" else value
    except Exception:
        return default


def infer_magnification(mpp: float) -> float | None:
    """Nearest standard objective for a pixel size (40x ~ 0.25 µm/px, 20x ~ 0.5 µm/px)."""
    if not mpp:
        return None
    return min((5, 10, 20, 40, 60, 80), key=lambda m: abs(math.log((10 / mpp) / m)))


def check_source(slide: Slide) -> tuple[CheckResult, SourceInfo]:
    """Check the largest file of the folder.

    Metadata problems (missing or unexpected magnification, pixel size, image type, …)
    are warnings: the file is still converted and the warning goes to the summary for
    review. Only problems that make conversion impossible are failures.
    """
    f = slide.largest
    res = CheckResult()
    info = SourceInfo(path=f.path)
    try:
        ds, info.pixel_data_offset = pixel_data_offset(f.path)
    except Exception as e:
        res.check(False, "Readable", f"cannot read DICOM header ({e})")
        return res, info

    sop = safe(lambda: ds.SOPClassUID)
    res.check(sop == WSI_SOP_CLASS, "Whole-slide image",
              "VL Whole Slide Microscopy" if sop == WSI_SOP_CLASS
              else f"SOP class is {getattr(sop, 'name', sop) or 'missing'}", fatal=False)

    image_type = [str(t).upper() for t in safe(lambda: list(ds.ImageType), [])]
    full_res = (len(image_type) >= 3 and image_type[0] == "ORIGINAL" and image_type[2] == "VOLUME"
                and (len(image_type) < 4 or image_type[3] == "NONE"))
    res.check(full_res, "Full resolution",
              f"ImageType {'/'.join(image_type) or 'missing'}"
              + ("" if full_res else " (expected ORIGINAL/.../VOLUME/NONE)"), fatal=False)

    info.width = int(safe(lambda: ds.TotalPixelMatrixColumns, 0))
    info.height = int(safe(lambda: ds.TotalPixelMatrixRows, 0))
    others = [s for s in slide.files if s is not f]
    bigger = [s for s in others if s.width * s.height > info.width * info.height]
    res.check(not bigger, "Largest level",
              f"{info.width:,} x {info.height:,} px is the largest of {len(others) + 1} "
              f"DICOM file(s) in the folder" if not bigger
              else f"{bigger[0].path.name} has more pixels ({bigger[0].width:,} x "
                   f"{bigger[0].height:,}) but a smaller file size", fatal=False)

    info.objective = safe(lambda: float(ds.OpticalPathSequence[0].ObjectiveLensPower))
    res.check(info.objective == EXPECTED_OBJECTIVE, "Objective 40x",
              f"objective lens power {info.objective:g}x" if info.objective is not None
              else "objective lens power not recorded", fatal=False)

    spacing = safe(lambda: [float(v) * 1000 for v in  # mm -> µm; (row, column) order
                            ds.SharedFunctionalGroupsSequence[0].PixelMeasuresSequence[0].PixelSpacing])
    if spacing is None:
        spacing = safe(lambda: [float(v) * 1000 for v in ds.PixelSpacing])
    mpp_y, info.mpp = (spacing if spacing and len(spacing) == 2 and min(spacing) > 0 else (0.0, 0.0))
    lo, hi = MPP_RANGE_40X
    res.check(lo <= info.mpp <= hi and abs(info.mpp - mpp_y) < 1e-3, "40x pixel size",
              f"{info.mpp:.4f} µm/px (40x range {lo}-{hi})"
              + (f", but {mpp_y:.4f} µm/px vertically" if info.mpp and abs(info.mpp - mpp_y) >= 1e-3 else "")
              if info.mpp else "pixel spacing not recorded", fatal=False)

    # Magnification written to the .svs: the recorded objective, else inferred from MPP.
    if info.objective is not None:
        info.app_mag, info.mag_source = info.objective, "metadata"
    elif info.mpp:
        info.app_mag, info.mag_source = infer_magnification(info.mpp), "inferred from MPP"
        res.warnings.append(f"Magnification: .svs labelled {info.app_mag:g}x, inferred from "
                            f"{info.mpp:.4f} µm/px")
    else:
        res.warnings.append("Magnification: .svs has no magnification or MPP; "
                            "set them in the viewer before measuring")
    if info.objective is not None and info.mpp:
        expected = infer_magnification(info.mpp)
        if expected != info.objective:
            res.warnings.append(f"Magnification: objective says {info.objective:g}x but "
                                f"{info.mpp:.4f} µm/px looks like {expected:g}x")

    # Everything below is needed to copy the tiles, so problems here stop the conversion.
    syntax = safe(lambda: ds.file_meta.TransferSyntaxUID)
    info.tile_w, info.tile_h = int(safe(lambda: ds.Columns, 0)), int(safe(lambda: ds.Rows, 0))
    info.frames = int(safe(lambda: ds.NumberOfFrames, 1))
    res.check(info.width > 0 and info.height > 0, "Image size",
              f"{info.width:,} x {info.height:,} px" if info.width and info.height
              else "TotalPixelMatrixColumns/Rows missing, so the tiles cannot be placed")
    expected_frames = (math.ceil(info.width / max(info.tile_w, 1))
                       * math.ceil(info.height / max(info.tile_h, 1)))
    layout_ok = (
        syntax == JPEG_BASELINE
        and safe(lambda: ds.DimensionOrganizationType, "TILED_FULL") == "TILED_FULL"
        and int(safe(lambda: ds.TotalPixelMatrixFocalPlanes, 1)) == 1
        and int(safe(lambda: ds.NumberOfOpticalPaths, 1)) == 1
        and int(safe(lambda: ds.SamplesPerPixel, 0)) == 3
        and int(safe(lambda: ds.BitsAllocated, 0)) == 8
        and info.frames == expected_frames > 0
    )
    if layout_ok:
        info.jpeg_colorspace, info.subsampling = safe(
            lambda: jpeg_colour_info(first_frame(f.path, info.pixel_data_offset)), ("unknown", None))
        layout_ok = info.jpeg_colorspace != "unknown"
    res.check(layout_ok, "Tile encoding",
              f"{info.frames:,} JPEG tiles of {info.tile_w}x{info.tile_h}, single focal plane"
              if layout_ok else
              f"unsupported: {getattr(syntax, 'name', syntax)}, "
              f"{safe(lambda: ds.DimensionOrganizationType, '?')}, {info.frames} frames "
              f"(this tool copies TILED_FULL JPEG-baseline tiles)")

    info.manufacturer = str(safe(lambda: ds.Manufacturer, ""))
    info.model = str(safe(lambda: ds.ManufacturerModelName, ""))
    info.icc = safe(lambda: bytes(ds.OpticalPathSequence[0].ICCProfile))
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
        *([f"AppMag = {info.app_mag:g}"] if info.app_mag else []),
        *([f"MPP = {info.mpp:.6f}"] if info.mpp else []),
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
        power = p.get("openslide.objective-power")
        if info.app_mag:
            res.check(power == f"{info.app_mag:g}", "Magnification kept",
                      f"objective power {power} ({info.mag_source})")
        mpp = float(p.get("openslide.mpp-x", 0) or 0)
        if info.mpp:
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
    notes: list[str] = field(default_factory=list)     # errors and reasons for skipping
    warnings: list[str] = field(default_factory=list)  # metadata issues to review


def needs_review(o: Outcome) -> bool:
    return bool(o.warnings) or o.status in ("skipped", "failed")


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
    for line in res.warnings:
        tqdm.write(f"    [WARN] {line}")
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
            (f"{i.app_mag:g}x" + ("*" if i.mag_source != "metadata" else "")) if i and i.app_mag else "-",
            f"{i.mpp:.4f}" if i and i.mpp else "-",
            human(o.slide.largest.size),
            human(o.output.stat().st_size) if o.output and o.output.exists() else "-",
            f"{o.seconds:.0f}s" if o.seconds else "-",
            f"{len(o.warnings)} warning(s)" if o.warnings else "-",
        ])
    header = ["Slide", "Status", "Pixels", "Mag", "MPP", "DICOM", "SVS", "Time", "Review"]
    widths = [max(len(str(r[k])) for r in rows + [header]) for k in range(len(header))]
    fmt = "  ".join(f"{{:<{w}}}" for w in widths)
    print(fmt.format(*header))
    print(fmt.format(*("-" * w for w in widths)))
    for r in rows:
        print(fmt.format(*r))

    if any(i and i.app_mag and i.mag_source != "metadata" for i in (o.info for o in outcomes)):
        print("* magnification inferred from pixel size (objective not recorded)")

    problems = [o for o in outcomes if o.notes or o.warnings]
    if problems:
        print("\nDetails:")
        for o in problems:
            for n in o.notes:
                print(f"  {o.slide.name}: {n}")
            for n in o.warnings:
                print(f"  {o.slide.name}: WARNING {n}")

    counts: dict[str, int] = {}
    for o in outcomes:
        counts[o.status] = counts.get(o.status, 0) + 1
    review = sum(1 for o in outcomes if needs_review(o))
    print("\n" + ", ".join(f"{v} {k}" for k, v in counts.items())
          + f" — {len(outcomes)} slide(s) in {elapsed:.0f}s"
          + (f"; {review} need review (see the needs_review column)" if review else ""))

    if not target.exists():
        return None
    report = target / f"dicom2svs_summary_{datetime.now():%Y%m%d_%H%M%S}.csv"
    with open(report, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["slide", "status", "needs_review", "warning_count", "warnings", "errors",
                    "source_folder", "source_file", "dicom_files_in_folder", "output_file",
                    "width", "height", "objective_recorded", "svs_magnification",
                    "magnification_source", "mpp", "dicom_bytes", "svs_bytes", "seconds"])
        for o in outcomes:
            i = o.info
            out_ok = o.output and o.output.exists()
            w.writerow([
                o.slide.name, o.status, "yes" if needs_review(o) else "no", len(o.warnings),
                " | ".join(o.warnings), " | ".join(o.notes),
                o.slide.folder, o.slide.largest.path, len(o.slide.files),
                o.output if out_ok else "",
                i.width or "" if i else "", i.height or "" if i else "",
                f"{i.objective:g}" if i and i.objective is not None else "",
                f"{i.app_mag:g}" if i and i.app_mag else "", i.mag_source if i else "",
                f"{i.mpp:.6f}" if i and i.mpp else "",
                o.slide.largest.size, o.output.stat().st_size if out_ok else "",
                f"{o.seconds:.1f}"])
    print(f"Summary saved to {report}")
    return report


# --------------------------------------------------------------------------- #
# Main
# --------------------------------------------------------------------------- #

def main() -> int:
    ap = argparse.ArgumentParser(
        description="Convert the largest DICOM whole-slide image in each folder to a "
                    "pyramidal Aperio .svs file. Missing or unexpected metadata (e.g. no "
                    "objective power) is recorded as a warning in the summary, not skipped.")
    ap.add_argument("dirs", nargs="+", type=Path, help="directories to scan (recursively)")
    ap.add_argument("-o", "--output", type=Path, required=True, help="target directory for .svs files")
    ap.add_argument("--dry-run", action="store_true", help="scan and check only, do not convert")
    ap.add_argument("--overwrite", action="store_true", help="replace existing .svs files")
    ap.add_argument("--strict", action="store_true",
                    help="skip slides with metadata warnings instead of converting them")
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
    print(f"Found {len(slides)} folder(s) with {sum(len(s.files) for s in slides)} DICOM "
          f"file(s); converting the largest file of each folder")
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
        o = Outcome(slide, "", info, out_path, warnings=list(checks.warnings))
        outcomes.append(o)

        if not checks.ok:
            o.status, o.output = "skipped", None
            o.notes += [f"check failed - {x}" for x in checks.failed]
            continue
        if args.strict and checks.warnings:
            o.status, o.output = "skipped", None
            o.notes.append("skipped by --strict because of metadata warnings")
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
