# Core API Reference

`core/` has no QGIS dependency and no I/O beyond raster/vector reads — it is
usable standalone from a notebook or test, and is what the `processing/`
algorithms wrap. Runtime dependencies are numpy, pyproj, and rasterio only
(see `requirements.txt`); inside QGIS, numpy and GDAL already ship with the
application.

## `core.geometry`

AEQD projection and ray sampling.

- `aeqd_proj4(lon0, lat0)` — PROJ string for an AEQD projection centered on
  `(lon0, lat0)`.
- `transformer_to_aeqd(lon0, lat0)` / `transformer_from_aeqd(lon0, lat0)` —
  `pyproj.Transformer` to/from the local AEQD plane.
- `ray_points(azimuth_deg, max_range_m, step_m)` → `(s, x, y)` — range and
  local-plane coordinates along one ray. `s` starts at `step_m`, not 0 — the
  radar's own cell is not part of the propagation path.
- `native_step(dem_transform, dem_crs, lat0=None)` — native ground sampling
  distance of a DEM in meters. Requires `lat0` when `dem_crs` is geographic
  (degrees), since meters-per-degree-longitude varies with latitude.

## `core.dem`

DEM access with explicit nodata handling.

- `Dem(path, vertical_datum=None)` — opens a single-band elevation raster
  (via rasterio), reads the whole band into memory once (a sweep samples
  hundreds of azimuths per site; re-reading from disk each time was the
  original bottleneck).
    - `.sample_lonlat(lons, lats, nodata_policy="raise")` — sample elevation
      at WGS84 points, reprojecting internally. `nodata_policy`:
      `"raise"` (default, raises `DemNodataError` on any nodata sample) or
      `"block"` (returns `+inf` for nodata cells — terrain that fully blocks
      the beam). Never silently returns 0 or NaN.
    - `.transform`, `.crs`, `.close()`, context-manager support.
- `DemNodataError(ValueError)` — raised per the policy above.
- `raster_cell_centers_lonlat(transform, crs, height, width)` — WGS84
  lon/lat of every cell center of a raster grid; used by Discover/Network to
  turn a mask or importance raster's own grid into query points.

## `core.horizon`

Accumulated terrain-horizon angle and minimum visible height. See
[Methodology](methodology.md#horizon-angle-and-minimum-visible-height) for
the derivation.

- `apparent_terrain_height(z, s, h0, k, earth_radius=EARTH_RADIUS_M)`
- `horizon_angle(z, s, h0, k, earth_radius=EARTH_RADIUS_M)` — running-max
  `θ(s)`.
- `min_visible_height(z, s, h0, k, theta=None, earth_radius=EARTH_RADIUS_M)`
  — `H_min(s)`; reuses `theta` if already computed.

## `core.blockage`

Partial beam blockage (Bech et al. 2003). See
[Methodology](methodology.md#partial-beam-blockage-bech-et-al-2003).

- `beam_height(s, h0, elevation_rad, k, earth_radius=EARTH_RADIUS_M)`
- `half_power_radius(s, beamwidth_rad)` — `a(s)`.
- `beam_block_frac(y, a)` — per-bin `PBB`.
- `cum_beam_block_frac(pbb, axis=-1)` — `CBB`, running max along range.
  **Never publish `PBB` directly** — see module docstring.

## `core.engine`

The per-site sweep, factored out for reuse by Evaluate and the Fase 0
notebook.

- `sweep_site(dem, lon0, lat0, h0, k, elevation_deg, beamwidth_deg, max_range_m, azimuth_step_deg=1.0, range_step_m=None, nodata_policy="block")`
  — full 360° sweep. Returns a dict: `azimuths_deg`, `s_m` (shared range
  axis), `theta_rad`, `h_min_m`, `cbb`, `lons_deg`, `lats_deg` — each
  `(n_azimuths, n_ranges)` except `azimuths_deg`/`s_m`.
- `h_min_at_point(dem, radar_lon, radar_lat, h0, k, query_lon, query_lat, range_step_m=None, nodata_policy="block")`
  — exact `H_min` at one arbitrary point, along the exact bearing (used by
  Network instead of interpolating between `sweep_site`'s fixed azimuths).
  Returns `h0` if the query point is closer than one range step.
- `site_merit(dem, importance_band, importance_transform, importance_crs, importance_nodata, lon0, lat0, h0, k, elevation_deg, beamwidth_deg, max_range_m, azimuth_step_deg=1.0, thresholds=(1000.0, 3000.0), blocked_sector_threshold=0.5)`
  — importance-weighted merit + largest blocked sector for one candidate;
  the two headline numbers the Robustness matrix sweeps over. Raises
  `ValueError` if the candidate's sweep doesn't overlap the importance layer
  at all — callers sweeping many candidates should catch this per candidate.

## `core.importance`

- `sample_importance_along_rays(band, transform, crs, nodata, lons, lats)` —
  importance value at each `(lon, lat)`; `0.0` outside the raster extent or
  at nodata (unlike `core.dem`, missing importance is a modeling choice, not
  a silent-bug risk, so it doesn't raise).

## `core.metrics`

Merit figures and penalties — plain array reductions, no I/O. Callers supply
`h_min`/`importance` on a shared grid and an explicit `mask` where needed.

- `weighted_visible_fraction(h_min, importance, threshold_m, mask=None)`
- `merit_scores(h_min, importance, thresholds=(1000.0, 3000.0), mask=None)`
  — dict of `{threshold: weighted_visible_fraction}`, deliberately never
  collapsed into one score.
- `largest_contiguous_blocked_sector_deg(cbb_by_azimuth, azimuth_step_deg, threshold=0.5)`
  — azimuths treated as circular.
- `network_residual_gap_fraction(h_min_stack, importance, threshold_m=3000.0, mask=None, axis=0)`
  — area seen by no radar in a stacked (`n_radars`, ...) `h_min` array.

All raise `ValueError` (via the private `_weighted_fraction`) if the total
importance weight over valid cells is zero or non-finite — there is no
default that silently returns a meaningless `0/0`.

## `core.manifest`

Reproducibility manifest. See
[Methodology](methodology.md#reproducibility-manifest).

- `sha256_file(path, chunk_size=1<<20)` — streamed digest, safe for
  multi-GB DEMs.
- `plugin_version(metadata_path=...)` — reads `metadata.txt`'s `version=`
  line.
- `gdal_version()` — via `rasterio.__gdal_version__`; inside QGIS this is
  the application's own bundled GDAL.
- `build_manifest(params, input_layers, plugin_version_override=None)` —
  assembles `{timestamp_utc, plugin_version, numpy_version, gdal_version,
  parameters, input_layers}`, hashing every path in `input_layers`.
- `write_manifest(path, manifest)` — pretty-printed JSON.
