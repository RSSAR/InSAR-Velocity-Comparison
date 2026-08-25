from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import rasterio
from rasterio.transform import from_origin


SCRIPT = Path(__file__).resolve().parents[1] / "compare_insar_velocity.py"


def _write_raster(path: Path, data: np.ndarray, transform, crs="EPSG:32611") -> None:
    profile = {
        "driver": "GTiff",
        "height": data.shape[0],
        "width": data.shape[1],
        "count": 1,
        "dtype": "float32",
        "crs": crs,
        "transform": transform,
        "nodata": np.nan,
    }
    with rasterio.open(path, "w", **profile) as dst:
        dst.write(data.astype("float32"), 1)


def test_median_reference_recovers_known_offset(tmp_path: Path) -> None:
    h, w = 60, 70
    transform = from_origin(500000, 3800000, 30, 30)
    rows, cols = np.indices((h, w))
    a = 0.05 * cols - 0.03 * rows
    b = a + 6.25

    a_path = tmp_path / "a.tif"
    b_path = tmp_path / "b.tif"
    outdir = tmp_path / "out"
    _write_raster(a_path, a, transform)
    _write_raster(b_path, b, transform)

    result = subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            str(a_path),
            str(b_path),
            "-o",
            str(outdir),
            "--reference",
            "median",
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    assert "Comparison complete" in result.stdout

    metrics = json.loads((outdir / "metrics.json").read_text())
    assert abs(metrics["reference_offset_subtracted_from_b"] - 6.25) < 1e-5
    assert abs(metrics["metrics_after_reference"]["bias_b_minus_a"]) < 1e-5
    assert metrics["metrics_after_reference"]["rmse"] < 1e-5


def test_stable_raster_reference(tmp_path: Path) -> None:
    h, w = 50, 50
    transform = from_origin(400000, 4100000, 20, 20)
    rows, cols = np.indices((h, w))
    a = 0.02 * cols + 0.01 * rows
    b = a + 4.0

    stable = np.zeros((h, w), dtype=np.float32)
    stable[10:30, 15:35] = 1.0

    a_path = tmp_path / "a.tif"
    b_path = tmp_path / "b.tif"
    stable_path = tmp_path / "stable.tif"
    outdir = tmp_path / "out_stable"
    _write_raster(a_path, a, transform)
    _write_raster(b_path, b, transform)
    _write_raster(stable_path, stable, transform)

    subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            str(a_path),
            str(b_path),
            "-o",
            str(outdir),
            "--reference",
            "stable",
            "--stable-area",
            str(stable_path),
            "--min-reference-pixels",
            "20",
        ],
        check=True,
        capture_output=True,
        text=True,
    )

    metrics = json.loads((outdir / "metrics.json").read_text())
    assert abs(metrics["reference_offset_subtracted_from_b"] - 4.0) < 1e-5
    assert abs(metrics["metrics_after_reference"]["median_bias_b_minus_a"]) < 1e-5
