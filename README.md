# InSAR Velocity Comparison

A practical command-line tool for **quantitatively comparing two geocoded InSAR velocity rasters** from different sensors, processing strategies, or software chains.

Typical use cases include:

- NISAR vs. Sentinel-1 velocity comparison;
- SBAS vs. DS-InSAR;
- different phase-linking methods;
- different atmospheric corrections;
- different reference-point or reference-area choices;
- reprocessed vs. previous InSAR solutions.

The script automatically finds the spatial overlap, aligns the two rasters onto a common grid, keeps common valid pixels, optionally harmonizes the reference frame, computes agreement statistics, and writes georeferenced difference products.

## Main outputs

The default difference convention is:

```text
Difference = B_corrected - A
```

When reference harmonization is enabled:

```text
B_corrected = B - reference_offset
```

The tool writes:

```text
difference_B_minus_A.tif
common_valid_mask.tif
metrics.json
metrics.csv
report.txt
```

Optional outputs include aligned copies of both rasters, the reference mask, and a PNG quicklook.

## Features

- Supports georeferenced rasters readable by Rasterio, especially GeoTIFF.
- Automatically detects the spatial intersection between the two velocity maps.
- Handles different:
  - raster extents;
  - pixel sizes;
  - raster grids;
  - coordinate reference systems.
- Reprojects one raster onto the selected target grid.
- Supports `nearest`, `bilinear`, and `cubic` resampling.
- Supports different bands for multiband rasters.
- Supports scale and offset conversion for each raster.
- Supports optional value-range filtering.
- Supports three reference modes:
  - `none` — compare original values;
  - `median` — remove the median B-A offset over all common valid pixels;
  - `stable` — remove the median B-A offset over a user-defined stable area.
- Stable areas can be supplied as:
  - polygon vector data (`.shp`, `.gpkg`, `.geojson`);
  - nonzero raster masks.
- Computes:
  - Bias;
  - median bias;
  - MAE;
  - RMSE;
  - standard deviation of the difference;
  - NMAD;
  - Pearson correlation coefficient;
  - linear-regression slope/intercept;
  - mean and median velocity for both rasters.
- Writes quantitative results before and after reference harmonization.
- Optional quicklook containing both velocity maps, the difference map, and a scatter plot.

## Requirements

- Python 3.9+
- NumPy
- Rasterio
- Fiona (only needed when `--stable-area` is a vector file)
- Matplotlib (only needed for `--plot`)

Install the required packages with:

```bash
pip install -r requirements.txt
```

## Basic usage

```bash
python compare_insar_velocity.py velocity_A.tif velocity_B.tif
```

Results are written to:

```text
velocity_comparison/
```

## Example: NISAR vs Sentinel-1

Assume both velocity products are already geocoded and expressed in `mm/yr`:

```bash
python compare_insar_velocity.py \
    sentinel1_velocity.tif \
    nisar_velocity.tif \
    -o S1_vs_NISAR \
    --reference median \
    --write-aligned \
    --plot
```

This aligns the NISAR velocity map to the Sentinel-1 grid by default and estimates a robust global reference offset using the median common difference.

## Stable-area reference harmonization

For a more defensible comparison, use a known stable region rather than the complete scene.

### Vector stable area

```bash
python compare_insar_velocity.py \
    sentinel1_velocity.tif \
    nisar_velocity.tif \
    -o S1_vs_NISAR \
    --reference stable \
    --stable-area stable_area.gpkg \
    --write-aligned \
    --plot
```

The reference offset is:

```text
median(B - A) within the stable area
```

and the corrected raster is:

```text
B_corrected = B - reference_offset
```

### Raster stable mask

A raster can also be used. Any finite nonzero pixel is considered part of the stable area:

```bash
python compare_insar_velocity.py A.tif B.tif \
    --reference stable \
    --stable-area stable_mask.tif
```

## Unit conversion

If one product is in `m/yr` and you want to compare in `mm/yr`:

```bash
python compare_insar_velocity.py A_m_per_yr.tif B_mm_per_yr.tif \
    --scale-a 1000 \
    --unit mm/yr
```

If LOS sign conventions are opposite, a scale factor of `-1` can also be used:

```bash
python compare_insar_velocity.py A.tif B.tif \
    --scale-b -1
```

This is useful when two processing chains use opposite positive-LOS conventions.

## Value filtering

To exclude clearly unrealistic or unwanted velocity values from both products:

```bash
python compare_insar_velocity.py A.tif B.tif \
    --valid-range -500 500
```

Only pixels where **both** products fall inside the inclusive range are retained.

## Target grid

By default, output products use raster A's grid:

```bash
--target-grid a
```

To use raster B's grid:

```bash
--target-grid b
```

## Resampling

Default:

```bash
--resampling bilinear
```

Available choices:

```text
nearest
bilinear
cubic
```

For continuous velocity fields, `bilinear` is generally appropriate. For discrete masks or categorical values, use `nearest`.

## Quantitative metrics

The following statistics are calculated on common valid pixels.

### Bias

```text
mean(B - A)
```

A positive bias means B is, on average, more positive than A.

### Mean absolute error (MAE)

```text
mean(|B - A|)
```

### Root mean square error (RMSE)

```text
sqrt(mean((B - A)^2))
```

### Difference standard deviation

Measures the spatial dispersion of the B-A residual field after any reference correction.

### NMAD

The normalized median absolute deviation is calculated as:

```text
1.4826 * median(|d - median(d)|)
```

where `d = B - A`.

NMAD is less sensitive to extreme outliers than standard deviation.

### Pearson correlation

Measures the linear spatial correspondence between the two velocity maps.

### Linear regression

The script also reports the least-squares relationship:

```text
B = slope * A + intercept
```

## Output files

### `difference_B_minus_A.tif`

Georeferenced residual velocity map:

```text
B_corrected - A
```

### `common_valid_mask.tif`

Binary mask of pixels used in the comparison.

### `reference_mask.tif`

Created when `--reference median` or `--reference stable` is used.

### `aligned_A.tif`

Created with `--write-aligned`.

### `aligned_B_corrected.tif`

Created with `--write-aligned`. This is B after spatial alignment and reference-offset correction.

### `metrics.json`

Machine-readable metadata and statistics.

### `metrics.csv`

Compact table containing metrics before and after reference correction.

### `report.txt`

Human-readable summary containing input paths, CRS, resolution, overlap bounds, reference offset, and comparison statistics.

### `comparison_quicklook.png`

Created with `--plot` and contains:

- Raster A;
- corrected Raster B;
- B-A difference;
- scatter plot with a 1:1 reference line.

## Recommended workflow for cross-sensor comparison

For NISAR and Sentinel-1, a defensible comparison workflow is:

1. Geocode both velocity products.
2. Convert both to the same physical unit.
3. Confirm that LOS sign conventions are consistent.
4. Use the same or physically comparable spatial region.
5. Define an independent stable area if possible.
6. Run this tool with `--reference stable`.
7. Inspect the residual GeoTIFF rather than relying only on a single correlation coefficient.
8. Report Bias/MAE/RMSE together with correlation and the number of common pixels.

A high correlation alone does not imply that the two velocity fields have the same absolute reference or magnitude.

## Complete example

```bash
python compare_insar_velocity.py \
    Sentinel1_velocity.tif \
    NISAR_velocity.tif \
    -o comparison_S1_NISAR \
    --target-grid a \
    --resampling bilinear \
    --unit mm/yr \
    --valid-range -300 300 \
    --reference stable \
    --stable-area stable_area.gpkg \
    --min-reference-pixels 100 \
    --write-aligned \
    --plot
```

## Important interpretation notes

This script performs **spatial and reference-frame harmonization**, but it does not make two sensors physically identical.

Differences may still arise from:

- different radar wavelengths;
- different incidence/heading geometry;
- different LOS projection vectors;
- different acquisition epochs and temporal sampling;
- atmospheric residuals;
- unwrapping errors;
- phase-linking strategy;
- spatial filtering or multilooking;
- reference-point/reference-area selection;
- different sensitivity to vegetation and decorrelation;
- real temporal changes in deformation rate.

Therefore, direct C-band vs L-band comparison should ideally be performed after accounting for geometry and temporal sampling where necessary.

## Command-line help

```bash
python compare_insar_velocity.py -h
```

## Version

Current version: `1.0.0`

## Author

Shuai Wang  
China University of Mining and Technology

## License

MIT License. See `LICENSE`.
