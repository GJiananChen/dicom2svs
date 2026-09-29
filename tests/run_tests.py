# /// script
# requires-python = ">=3.10"
# dependencies = []
# ///
"""
End-to-end test for dicom2svs on synthetic slides. Works on Windows, Linux and macOS.

  uv run tests/run_tests.py
"""

import csv
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent

EXPECTED = {
    # slide: (status, needs_review)
    "good_40x": ("converted", "no"),
    "no_metadata": ("converted", "yes"),
    "jpeg2000": ("skipped", "yes"),
}


def run(*args: str) -> subprocess.CompletedProcess:
    print("$", " ".join(args), flush=True)
    return subprocess.run(args, text=True, encoding="utf-8", errors="replace",
                          stdout=subprocess.PIPE, stderr=subprocess.STDOUT)


def main() -> int:
    work = Path(tempfile.mkdtemp(prefix="dicom2svs_test_"))
    try:
        slides, out = work / "slides", work / "svs out"  # space in path on purpose
        r = run("uv", "run", str(REPO / "tests" / "make_test_slides.py"), str(slides))
        print(r.stdout)
        if r.returncode:
            return r.returncode

        r = run("uv", "run", str(REPO / "dicom2svs.py"), str(slides), "-o", str(out))
        print(r.stdout[-6000:])
        failures = []
        if r.returncode != 0:
            failures.append(f"dicom2svs exited with {r.returncode}")

        reports = sorted(out.glob("dicom2svs_summary_*.csv"))
        if not reports:
            failures.append("no summary CSV written")
        else:
            with open(reports[-1], encoding="utf-8-sig", newline="") as fh:
                rows = {row["slide"]: row for row in csv.DictReader(fh)}
            for slide, (status, review) in EXPECTED.items():
                row = rows.get(slide)
                if row is None:
                    failures.append(f"{slide}: missing from summary")
                    continue
                got = (row["status"], row["needs_review"])
                if got != (status, review):
                    failures.append(f"{slide}: expected {(status, review)}, got {got} "
                                    f"(errors: {row['errors']})")
            for slide in ("good_40x", "no_metadata"):
                if not (out / f"{slide}.svs").is_file():
                    failures.append(f"{slide}.svs not written")
            if (out / "jpeg2000.svs").exists():
                failures.append("jpeg2000.svs should not have been written")
            if rows.get("good_40x", {}).get("svs_magnification") != "40":
                failures.append("good_40x.svs should be labelled 40x")

        # Second run must leave existing files alone.
        r = run("uv", "run", str(REPO / "dicom2svs.py"), str(slides), "-o", str(out))
        if "2 exists" not in r.stdout:
            failures.append("rerun did not report existing outputs")

        if failures:
            print("\nFAILED:\n  " + "\n  ".join(failures))
            return 1
        print("\nAll tests passed.")
        return 0
    finally:
        shutil.rmtree(work, ignore_errors=True)


if __name__ == "__main__":
    sys.exit(main())
