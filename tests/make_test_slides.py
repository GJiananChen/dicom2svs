# /// script
# requires-python = ">=3.10"
# dependencies = ["pydicom>=3.0", "pillow>=10", "numpy"]
# ///
"""
Write small synthetic DICOM whole-slide images for testing dicom2svs (no patient data).

  uv run tests/make_test_slides.py OUT_DIR

Creates:
  OUT_DIR/good_40x/        full-resolution 40x slide + a 1/4 resampled level
  OUT_DIR/no_metadata/     same slide with objective, pixel spacing and ImageType removed
  OUT_DIR/jpeg2000/        a slide with an unsupported transfer syntax (must be skipped)
"""

import io
import sys
from pathlib import Path

import numpy as np
import pydicom
from PIL import Image
from pydicom.dataset import Dataset, FileMetaDataset
from pydicom.encaps import encapsulate
from pydicom.uid import JPEGBaseline8Bit, JPEG2000Lossless, generate_uid

WSI_SOP_CLASS = "1.2.840.10008.5.1.4.1.1.77.1.6"
TILE = 256


def synthetic_image(width: int, height: int) -> np.ndarray:
    """A pink/purple H&E-like pattern with sharp detail, so tile mix-ups are detectable."""
    y, x = np.mgrid[0:height, 0:width]
    r = 200 + 40 * np.sin(x / 37.0)
    g = 120 + 60 * np.cos(y / 23.0)
    b = 180 + 50 * np.sin((x + y) / 51.0)
    img = np.stack([r, g, b], axis=-1)
    img[(x // 64 + y // 64) % 7 == 0] = (90, 40, 140)  # "nuclei" blocks
    return np.clip(img, 0, 255).astype(np.uint8)


def tiles(img: np.ndarray) -> list[bytes]:
    h, w, _ = img.shape
    frames = []
    for ty in range(0, h, TILE):
        for tx in range(0, w, TILE):
            tile = np.full((TILE, TILE, 3), 255, np.uint8)
            part = img[ty:ty + TILE, tx:tx + TILE]
            tile[:part.shape[0], :part.shape[1]] = part
            buf = io.BytesIO()
            Image.fromarray(tile).save(buf, "JPEG", quality=90, subsampling=2)
            frames.append(buf.getvalue())
    return frames


def write_level(path: Path, img: np.ndarray, series_uid: str, mpp: float, image_type: list[str],
                objective: float | None = 40, syntax=JPEGBaseline8Bit) -> None:
    h, w, _ = img.shape
    frames = tiles(img)

    meta = FileMetaDataset()
    meta.MediaStorageSOPClassUID = WSI_SOP_CLASS
    meta.MediaStorageSOPInstanceUID = generate_uid()
    meta.TransferSyntaxUID = syntax

    ds = Dataset()
    ds.file_meta = meta
    ds.SOPClassUID = WSI_SOP_CLASS
    ds.SOPInstanceUID = meta.MediaStorageSOPInstanceUID
    ds.StudyInstanceUID = generate_uid()
    ds.SeriesInstanceUID = series_uid
    ds.Modality = "SM"
    ds.Manufacturer = "Synthetic"
    ds.ManufacturerModelName = "dicom2svs test"
    ds.PatientName = "TEST^SLIDE"
    ds.PatientID = "TEST0001"
    if image_type:
        ds.ImageType = image_type
    ds.SamplesPerPixel = 3
    ds.PhotometricInterpretation = "YBR_FULL_422"
    ds.PlanarConfiguration = 0
    ds.BitsAllocated = ds.BitsStored = 8
    ds.HighBit = 7
    ds.PixelRepresentation = 0
    ds.Rows = ds.Columns = TILE
    ds.NumberOfFrames = len(frames)
    ds.TotalPixelMatrixColumns = w
    ds.TotalPixelMatrixRows = h
    ds.TotalPixelMatrixFocalPlanes = 1
    ds.NumberOfOpticalPaths = 1
    ds.DimensionOrganizationType = "TILED_FULL"

    path_item = Dataset()
    path_item.OpticalPathIdentifier = "0"
    if objective is not None:
        path_item.ObjectiveLensPower = objective
    ds.OpticalPathSequence = [path_item]

    shared = Dataset()
    if mpp:
        measures = Dataset()
        measures.PixelSpacing = [mpp / 1000, mpp / 1000]  # mm
        measures.SliceThickness = 0.0017
        shared.PixelMeasuresSequence = [measures]
    ds.SharedFunctionalGroupsSequence = [shared]

    ds.PixelData = encapsulate(frames)
    ds["PixelData"].VR = "OB"
    ds["PixelData"].is_undefined_length = True
    path.parent.mkdir(parents=True, exist_ok=True)
    ds.save_as(path, enforce_file_format=True)


def main() -> None:
    out = Path(sys.argv[1] if len(sys.argv) > 1 else "test_slides")
    base = synthetic_image(5000, 3000)
    low = np.asarray(Image.fromarray(base).reduce(4))

    uid = generate_uid()
    write_level(out / "good_40x" / "slide_1_1.dcm", base, uid, 0.25,
                ["ORIGINAL", "PRIMARY", "VOLUME", "NONE"])
    write_level(out / "good_40x" / "slide_1_2.dcm", low, uid, 1.0,
                ["DERIVED", "PRIMARY", "VOLUME", "RESAMPLED"])

    write_level(out / "no_metadata" / "slide.dcm", base, generate_uid(), 0, [], objective=None)

    # Tiles are JPEG baseline but the header claims JPEG 2000: the tool must skip it.
    write_level(out / "jpeg2000" / "slide.dcm", base, generate_uid(), 0.25,
                ["ORIGINAL", "PRIMARY", "VOLUME", "NONE"], syntax=JPEG2000Lossless)
    print(f"Test slides written to {out.resolve()}")


if __name__ == "__main__":
    main()
