#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Compare two geocoded InSAR velocity rasters on a common grid.

The tool aligns the rasters spatially, keeps their common valid pixels,
optionally harmonizes their reference frames, calculates quantitative
agreement metrics, and writes georeferenced comparison products.

Difference convention:
    difference = B_corrected - A

Reference correction convention:
    B_corrected = B - reference_offset

where reference_offset is the median of (B - A) over either the complete
common valid area or a user-supplied stable area.

Author: Shuai Wang
Affiliation: China University of Mining and Technology
Version: 1.0.0
Release date: 2026-08-25
Repository: https://github.com/RSSAR/InSAR-Velocity-Comparison
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import rasterio
from rasterio.enums import Resampling
from rasterio.features import geometry_mask
from rasterio.windows import Window, from_bounds, transform as window_transform
from rasterio.warp import reproject, transform_bounds, transform_geom


VERSION = "1.0.0"


@dataclass
class GridData:
    a: np.ndarray
    b: np.ndarray
    transform: rasterio.Affine
    crs: Any
    profile: dict[str, Any]
    bounds: tuple[float, float, float, float]
    resolution: tuple[float, float]
    target_grid: str


def finite_mask(arr: np.ndarray) -> np.ndarray:
    return np.isfinite(arr)


def read_float_band(ds: rasterio.DatasetReader, band: int, window: Window | None = None) -> np.ndarray:
    data = ds.read(band, window=window, masked=True).astype(np.float64)
    return np.asarray(data.filled(np.nan), dtype=np.float64)


def _aligned_window(ds: rasterio.DatasetReader, bounds: tuple[float, float, float, float]) -> Window:
    win = from_bounds(*bounds, transform=ds.transform)
    col0 = max(0, math.floor(win.col_off))
    row0 = max(0, math.floor(win.row_off))
    col1 = min(ds.width, math.ceil(win.col_off + win.width))
    row1 = min(ds.height, math.ceil(win.row_off + win.height))
    if col1 <= col0 or row1 <= row0:
        raise ValueError("The rasters do not have a usable common overlap.")
    return Window(col0, row0, col1 - col0, row1 - row0)


def _intersect_bounds(
    ref_ds: rasterio.DatasetReader,
    other_ds: rasterio.DatasetReader,
) -> tuple[float, float, float, float]:
    if ref_ds.crs is None or other_ds.crs is None:
        raise ValueError("Both rasters must have a valid CRS.")

    other_in_ref = transform_bounds(
        other_ds.crs,
        ref_ds.crs,
        *other_ds.bounds,
        densify_pts=21,
    )
    left = max(ref_ds.bounds.left, other_in_ref[0])
    bottom = max(ref_ds.bounds.bottom, other_in_ref[1])
    right = min(ref_ds.bounds.right, other_in_ref[2])
    top = min(ref_ds.bounds.top, other_in_ref[3])
    if left >= right or bottom >= top:
        raise ValueError("The two rasters do not overlap spatially.")
    return left, bottom, right, top


def _resampling(name: str) -> Resampling:
    mapping = {
        "nearest": Resampling.nearest,
        "bilinear": Resampling.bilinear,
        "cubic": Resampling.cubic,
    }
    return mapping[name]


def prepare_common_grid(
    path_a: str,
    path_b: str,
    band_a: int,
    band_b: int,
    target_grid: str,
    resampling: str,
    scale_a: float,
    scale_b: float,
    offset_a: float,
    offset_b: float,
) -> GridData:
    with rasterio.open(path_a) as a_ds, rasterio.open(path_b) as b_ds:
        if band_a < 1 or band_a > a_ds.count:
            raise ValueError(f"--band-a must be between 1 and {a_ds.count}")
        if band_b < 1 or band_b > b_ds.count:
            raise ValueError(f"--band-b must be between 1 and {b_ds.count}")

        if target_grid == "a":
            ref_ds, other_ds = a_ds, b_ds
            ref_band, other_band = band_a, band_b
        else:
            ref_ds, other_ds = b_ds, a_ds
            ref_band, other_band = band_b, band_a

        overlap = _intersect_bounds(ref_ds, other_ds)
        win = _aligned_window(ref_ds, overlap)
        dst_transform = window_transform(win, ref_ds.transform)
        dst_height = int(win.height)
        dst_width = int(win.width)

        ref_arr = read_float_band(ref_ds, ref_band, win)
        other_src = read_float_band(other_ds, other_band)
        other_arr = np.full((dst_height, dst_width), np.nan, dtype=np.float64)

        reproject(
            source=other_src,
            destination=other_arr,
            src_transform=other_ds.transform,
            src_crs=other_ds.crs,
            src_nodata=np.nan,
            dst_transform=dst_transform,
            dst_crs=ref_ds.crs,
            dst_nodata=np.nan,
            resampling=_resampling(resampling),
            init_dest_nodata=True,
        )

        if target_grid == "a":
            arr_a = ref_arr * scale_a + offset_a
            arr_b = other_arr * scale_b + offset_b
        else:
            arr_b = ref_arr * scale_b + offset_b
            arr_a = other_arr * scale_a + offset_a

        profile = ref_ds.profile.copy()
        profile.update(
            driver="GTiff",
            count=1,
            height=dst_height,
            width=dst_width,
            transform=dst_transform,
            crs=ref_ds.crs,
            dtype="float32",
            nodata=np.nan,
            compress="deflate",
            predictor=3,
            BIGTIFF="IF_SAFER",
        )

        left, bottom, right, top = rasterio.transform.array_bounds(
            dst_height, dst_width, dst_transform
        )
        return GridData(
            a=arr_a,
            b=arr_b,
            transform=dst_transform,
            crs=ref_ds.crs,
            profile=profile,
            bounds=(left, bottom, right, top),
            resolution=(abs(dst_transform.a), abs(dst_transform.e)),
            target_grid=target_grid,
        )


def _read_stable_raster(
    path: str,
    grid: GridData,
) -> np.ndarray:
    with rasterio.open(path) as ds:
        src = read_float_band(ds, 1)
        dst = np.zeros(grid.a.shape, dtype=np.float32)
        reproject(
            source=np.where(np.isfinite(src) & (src != 0), 1.0, 0.0).astype(np.float32),
            destination=dst,
            src_transform=ds.transform,
            src_crs=ds.crs,
            dst_transform=grid.transform,
            dst_crs=grid.crs,
            src_nodata=0.0,
            dst_nodata=0.0,
            resampling=Resampling.nearest,
        )
    return dst > 0.5


def _read_stable_vector(
    path: str,
    grid: GridData,
) -> np.ndarray:
    try:
        import fiona
    except ImportError as exc:
        raise RuntimeError(
            "Vector stable areas require Fiona. Install it or provide a raster mask."
        ) from exc

    shapes = []
    with fiona.open(path) as src:
        src_crs = src.crs_wkt or src.crs
        if not src_crs:
            raise ValueError("Stable-area vector has no CRS.")
        for feat in src:
            geom = feat.get("geometry")
            if geom is None:
                continue
            shapes.append(transform_geom(src_crs, grid.crs, geom, precision=-1))

    if not shapes:
        raise ValueError("Stable-area vector contains no usable geometries.")

    return geometry_mask(
        shapes,
        out_shape=grid.a.shape,
        transform=grid.transform,
        invert=True,
        all_touched=False,
    )


def stable_area_mask(path: str, grid: GridData) -> np.ndarray:
    suffix = Path(path).suffix.lower()
    vector_suffixes = {".shp", ".gpkg", ".geojson", ".json", ".kml"}
    if suffix in vector_suffixes:
        return _read_stable_vector(path, grid)
    return _read_stable_raster(path, grid)


def regression(a: np.ndarray, b: np.ndarray) -> tuple[float, float]:
    if a.size < 2 or np.nanstd(a) == 0:
        return float("nan"), float("nan")
    slope, intercept = np.polyfit(a, b, 1)
    return float(slope), float(intercept)


def metrics(a: np.ndarray, b: np.ndarray, mask: np.ndarray) -> dict[str, float | int]:
    valid = mask & finite_mask(a) & finite_mask(b)
    av = a[valid]
    bv = b[valid]
    if av.size == 0:
        raise ValueError("No common valid pixels are available for statistics.")

    diff = bv - av
    bias = float(np.mean(diff))
    median_bias = float(np.median(diff))
    mae = float(np.mean(np.abs(diff)))
    rmse = float(np.sqrt(np.mean(diff**2)))
    std = float(np.std(diff))
    nmad = float(1.4826 * np.median(np.abs(diff - np.median(diff))))

    if av.size >= 2 and np.std(av) > 0 and np.std(bv) > 0:
        corr = float(np.corrcoef(av, bv)[0, 1])
    else:
        corr = float("nan")

    slope, intercept = regression(av, bv)
    return {
        "n_valid": int(av.size),
        "mean_a": float(np.mean(av)),
        "mean_b": float(np.mean(bv)),
        "median_a": float(np.median(av)),
        "median_b": float(np.median(bv)),
        "bias_b_minus_a": bias,
        "median_bias_b_minus_a": median_bias,
        "mae": mae,
        "rmse": rmse,
        "std_difference": std,
        "nmad_difference": nmad,
        "pearson_r": corr,
        "regression_slope_b_vs_a": slope,
        "regression_intercept_b_vs_a": intercept,
    }


def write_float_tif(path: Path, data: np.ndarray, grid: GridData) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    out = np.where(np.isfinite(data), data, np.nan).astype(np.float32)
    with rasterio.open(path, "w", **grid.profile) as dst:
        dst.write(out, 1)


def write_mask_tif(path: Path, mask: np.ndarray, grid: GridData) -> None:
    profile = grid.profile.copy()
    profile.update(dtype="uint8", nodata=0, predictor=1)
    with rasterio.open(path, "w", **profile) as dst:
        dst.write(mask.astype(np.uint8), 1)


def _fmt(value: float | int) -> str:
    if isinstance(value, (int, np.integer)):
        return str(int(value))
    if not np.isfinite(value):
        return "nan"
    return f"{float(value):.6f}"


def write_metrics_csv(path: Path, before: dict[str, Any], after: dict[str, Any]) -> None:
    keys = list(before.keys())
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["metric", "before_reference", "after_reference"])
        for key in keys:
            writer.writerow([key, before[key], after[key]])


def write_report(
    path: Path,
    args: argparse.Namespace,
    grid: GridData,
    offset: float,
    common_mask: np.ndarray,
    reference_mask: np.ndarray | None,
    before: dict[str, Any],
    after: dict[str, Any],
) -> None:
    lines = [
        "InSAR Velocity Comparison Report",
        "=" * 32,
        f"Raster A: {Path(args.raster_a).resolve()}",
        f"Raster B: {Path(args.raster_b).resolve()}",
        f"Difference convention: B_corrected - A",
        f"Target grid: {grid.target_grid.upper()}",
        f"CRS: {grid.crs}",
        f"Resolution: {grid.resolution[0]:.12g}, {grid.resolution[1]:.12g}",
        "Overlap bounds (left, bottom, right, top): "
        + ", ".join(f"{x:.12g}" for x in grid.bounds),
        f"Resampling: {args.resampling}",
        f"Scale A / B: {args.scale_a} / {args.scale_b}",
        f"Offset A / B: {args.offset_a} / {args.offset_b}",
        f"Unit label: {args.unit}",
        f"Reference mode: {args.reference}",
        f"Reference offset subtracted from B: {offset:.12g} {args.unit}",
        f"Common valid pixels: {int(common_mask.sum())}",
    ]
    if args.stable_area:
        lines.append(f"Stable area: {Path(args.stable_area).resolve()}")
    if reference_mask is not None:
        lines.append(f"Reference pixels: {int(reference_mask.sum())}")
    if args.valid_range is not None:
        lines.append(f"Valid range: {args.valid_range[0]} to {args.valid_range[1]} {args.unit}")

    lines.extend(["", "Metrics", "-------", f"{'Metric':34s} {'Before':>16s} {'After':>16s}"])
    for key in before:
        lines.append(f"{key:34s} {_fmt(before[key]):>16s} {_fmt(after[key]):>16s}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_plot(
    path: Path,
    a: np.ndarray,
    b: np.ndarray,
    diff: np.ndarray,
    mask: np.ndarray,
    unit: str,
    max_scatter: int,
) -> None:
    try:
        import matplotlib.pyplot as plt
    except ImportError as exc:
        raise RuntimeError("--plot requires matplotlib.") from exc

    valid = mask & finite_mask(a) & finite_mask(b)
    rows, cols = np.nonzero(valid)
    if rows.size == 0:
        return

    rng = np.random.default_rng(42)
    if rows.size > max_scatter:
        idx = rng.choice(rows.size, size=max_scatter, replace=False)
        rows_s = rows[idx]
        cols_s = cols[idx]
    else:
        rows_s, cols_s = rows, cols

    fig, axes = plt.subplots(2, 2, figsize=(12, 10), constrained_layout=True)
    im0 = axes[0, 0].imshow(np.where(valid, a, np.nan))
    axes[0, 0].set_title("Raster A")
    fig.colorbar(im0, ax=axes[0, 0], label=unit)

    im1 = axes[0, 1].imshow(np.where(valid, b, np.nan))
    axes[0, 1].set_title("Raster B (reference-corrected)")
    fig.colorbar(im1, ax=axes[0, 1], label=unit)

    im2 = axes[1, 0].imshow(np.where(valid, diff, np.nan))
    axes[1, 0].set_title("Difference: B - A")
    fig.colorbar(im2, ax=axes[1, 0], label=unit)

    av = a[rows_s, cols_s]
    bv = b[rows_s, cols_s]
    axes[1, 1].scatter(av, bv, s=2, alpha=0.25)
    lo = float(np.nanmin([np.nanmin(av), np.nanmin(bv)]))
    hi = float(np.nanmax([np.nanmax(av), np.nanmax(bv)]))
    axes[1, 1].plot([lo, hi], [lo, hi], "k--", linewidth=1)
    axes[1, 1].set_xlabel(f"Raster A ({unit})")
    axes[1, 1].set_ylabel(f"Raster B ({unit})")
    axes[1, 1].set_title(f"Scatter (n={av.size:,})")
    axes[1, 1].set_aspect("equal", adjustable="box")

    for ax in axes.flat[:3]:
        ax.set_xticks([])
        ax.set_yticks([])

    fig.savefig(path, dpi=180)
    plt.close(fig)


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=(
            "Align and quantitatively compare two geocoded InSAR velocity rasters. "
            "Outputs use the difference convention B_corrected - A."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("raster_a", help="Reference comparison raster A (GeoTIFF or rasterio-readable raster).")
    p.add_argument("raster_b", help="Comparison raster B (GeoTIFF or rasterio-readable raster).")
    p.add_argument("-o", "--outdir", default="velocity_comparison", help="Output directory.")
    p.add_argument("--band-a", type=int, default=1, help="Band number for raster A.")
    p.add_argument("--band-b", type=int, default=1, help="Band number for raster B.")
    p.add_argument("--target-grid", choices=["a", "b"], default="a", help="Grid used for the common output raster.")
    p.add_argument("--resampling", choices=["nearest", "bilinear", "cubic"], default="bilinear", help="Resampling used during reprojection.")
    p.add_argument("--scale-a", type=float, default=1.0, help="Multiply raster A values by this factor.")
    p.add_argument("--scale-b", type=float, default=1.0, help="Multiply raster B values by this factor.")
    p.add_argument("--offset-a", type=float, default=0.0, help="Add this constant to raster A after scaling.")
    p.add_argument("--offset-b", type=float, default=0.0, help="Add this constant to raster B after scaling.")
    p.add_argument("--unit", default="mm/yr", help="Unit label used in reports and plots.")
    p.add_argument(
        "--valid-range",
        type=float,
        nargs=2,
        metavar=("MIN", "MAX"),
        help="Keep pixels only when both raster values lie within this inclusive range.",
    )
    p.add_argument(
        "--reference",
        choices=["none", "median", "stable"],
        default="none",
        help=(
            "Reference-frame harmonization. 'median' removes the median B-A over all common pixels; "
            "'stable' removes the median B-A inside --stable-area."
        ),
    )
    p.add_argument(
        "--stable-area",
        help="Stable-area polygon (.shp/.gpkg/.geojson) or nonzero raster mask, required for --reference stable.",
    )
    p.add_argument("--min-reference-pixels", type=int, default=20, help="Minimum valid pixels required for reference estimation.")
    p.add_argument("--write-aligned", action="store_true", help="Also write aligned A and corrected B rasters.")
    p.add_argument("--plot", action="store_true", help="Write a PNG quicklook with maps and scatter plot.")
    p.add_argument("--max-scatter", type=int, default=100000, help="Maximum points drawn in the scatter plot.")
    p.add_argument("--version", action="version", version=f"%(prog)s {VERSION}")
    return p


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)

    if args.reference == "stable" and not args.stable_area:
        raise ValueError("--reference stable requires --stable-area.")
    if args.min_reference_pixels < 1:
        raise ValueError("--min-reference-pixels must be >= 1.")
    if args.max_scatter < 1:
        raise ValueError("--max-scatter must be >= 1.")
    if args.valid_range is not None and args.valid_range[0] >= args.valid_range[1]:
        raise ValueError("--valid-range requires MIN < MAX.")

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    grid = prepare_common_grid(
        args.raster_a,
        args.raster_b,
        args.band_a,
        args.band_b,
        args.target_grid,
        args.resampling,
        args.scale_a,
        args.scale_b,
        args.offset_a,
        args.offset_b,
    )

    common = finite_mask(grid.a) & finite_mask(grid.b)
    if args.valid_range is not None:
        vmin, vmax = args.valid_range
        common &= (grid.a >= vmin) & (grid.a <= vmax) & (grid.b >= vmin) & (grid.b <= vmax)

    if not np.any(common):
        raise ValueError("No common valid pixels remain after masking/filtering.")

    before = metrics(grid.a, grid.b, common)
    reference_mask: np.ndarray | None = None
    offset = 0.0

    if args.reference == "median":
        reference_mask = common.copy()
    elif args.reference == "stable":
        stable = stable_area_mask(args.stable_area, grid)
        reference_mask = common & stable

    if reference_mask is not None:
        nref = int(reference_mask.sum())
        if nref < args.min_reference_pixels:
            raise ValueError(
                f"Only {nref} valid reference pixels are available; "
                f"at least {args.min_reference_pixels} are required."
            )
        offset = float(np.median((grid.b - grid.a)[reference_mask]))

    b_corrected = grid.b - offset
    after = metrics(grid.a, b_corrected, common)
    difference = np.where(common, b_corrected - grid.a, np.nan)

    write_float_tif(outdir / "difference_B_minus_A.tif", difference, grid)
    write_mask_tif(outdir / "common_valid_mask.tif", common, grid)

    if reference_mask is not None:
        write_mask_tif(outdir / "reference_mask.tif", reference_mask, grid)

    if args.write_aligned:
        write_float_tif(outdir / "aligned_A.tif", np.where(common, grid.a, np.nan), grid)
        write_float_tif(outdir / "aligned_B_corrected.tif", np.where(common, b_corrected, np.nan), grid)

    payload = {
        "version": VERSION,
        "raster_a": str(Path(args.raster_a).resolve()),
        "raster_b": str(Path(args.raster_b).resolve()),
        "difference_convention": "B_corrected - A",
        "target_grid": grid.target_grid,
        "crs": str(grid.crs),
        "resolution": list(grid.resolution),
        "overlap_bounds": list(grid.bounds),
        "reference_mode": args.reference,
        "reference_offset_subtracted_from_b": offset,
        "unit": args.unit,
        "metrics_before_reference": before,
        "metrics_after_reference": after,
    }
    (outdir / "metrics.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
    write_metrics_csv(outdir / "metrics.csv", before, after)
    write_report(
        outdir / "report.txt",
        args,
        grid,
        offset,
        common,
        reference_mask,
        before,
        after,
    )

    if args.plot:
        write_plot(
            outdir / "comparison_quicklook.png",
            grid.a,
            b_corrected,
            difference,
            common,
            args.unit,
            args.max_scatter,
        )

    print("Comparison complete")
    print(f"  valid pixels     : {int(common.sum()):,}")
    print(f"  reference mode   : {args.reference}")
    print(f"  reference offset : {offset:.6f} {args.unit}")
    print(f"  bias (B-A)       : {after['bias_b_minus_a']:.6f} {args.unit}")
    print(f"  MAE              : {after['mae']:.6f} {args.unit}")
    print(f"  RMSE             : {after['rmse']:.6f} {args.unit}")
    print(f"  Pearson r        : {after['pearson_r']:.6f}")
    print(f"  output directory : {outdir.resolve()}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(2)
