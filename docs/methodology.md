# Methodology

The full derivation lives in `Plan_Tecnico_Radar_Siting_Toolkit.md` (repo
root, not distributed with the plugin). This page summarizes the geometry
actually implemented in `core/`.

## Geometry: AEQD projection

All geometry for a single radar's analysis happens in a local azimuthal
equidistant (AEQD) plane centered on that radar (`core/geometry.py`): x/y in
meters, distance from the radar is exact Euclidean distance, azimuth from the
radar is exact bearing. This avoids the anisotropic distortion UTM introduces
at 250 km range and across zone boundaries.

Each ray is sampled from `step_m` out to `max_range_m`, at the DEM's own
native ground sampling distance by default (`native_step`) — never a coarser
fixed step, or intermediate ridges are skipped entirely. Because GLO-30 and
FABDEM ship in EPSG:4326 (degrees, not meters), `native_step` converts pixel
size to a geodesic distance at the radar's latitude; meters-per-degree
shrinks with latitude, so this cannot be a fixed constant.

## Horizon angle and minimum visible height

For each range `s` along a ray, `core/horizon.py` computes the apparent
terrain height relative to the radar, net of earth-curvature drop:

```
Δz(s) = z(s) − h0 − s² / (2 k R)
```

where `h0` is antenna height (ground elevation + tower height), `k` is the
effective earth-radius factor (1.0 = standard, 4/3 = typical tropospheric
refraction, higher = superrefractive), and `R` is Earth's radius.

The horizon angle is the running maximum of `atan2(Δz(s), s)` — the minimum
clear elevation angle at each range, since anything below it is blocked by
nearer terrain:

```
θ(s) = running_max( atan2(Δz(s), s) )
```

Minimum visible height above local terrain follows directly:

```
H_min(s) = h0 + s·tan(θ(s)) + s² / (2 k R) − z(s)
```

This is one `O(N)` pass per ray, replacing a discrete elevation-angle loop:
everything downstream (merit scores, network coverage) derives from this one
array.

## Partial beam blockage (Bech et al. 2003)

`core/blockage.py` reimplements the circular-segment PBB geometry of Bech et
al. (2003) in ~50 lines of plain numpy, dropping the `wradlib` runtime
dependency (scipy + matplotlib + xarray + xradar, conda-only install path —
not installable from the QGIS Plugin Manager). `wradlib` is kept only as an
offline validation reference (see `notebooks/`).

Beam center height at range `s`, for scan elevation `elevation_deg`:

```
beam_height(s) = h0 + s·tan(elevation_rad) + s² / (2 k R)
```

Half-power (−3 dB) beam radius at range `s`, for beamwidth `beamwidth_deg`:

```
a(s) = s·tan(beamwidth_rad / 2)
```

Per-bin blockage fraction `PBB` is the circular-segment overlap between the
beam's circular cross-section and the terrain, as a function of
`y = terrain_height − beam_center_height` and `a`: `0` when the beam is fully
clear (`y ≤ −a`), `1` when fully blocked (`y ≥ a`), the closed-form
circular-segment fraction in between.

**The cumulative form is mandatory for any published map.** `CBB(s)` is the
running maximum of `PBB` along range — a beam once blocked stays blocked at
every greater range regardless of what the terrain does next. Publishing raw
per-bin `PBB` instead produces wrong results across the whole sector behind
an obstacle.

## Figures of merit

`core/metrics.py` turns `H_min` and `CBB` into decision-grade numbers,
always importance-weighted rather than collapsed into a single score:

- **Weighted visible fraction** — importance-weighted fraction of the area
  where `H_min ≤ threshold`, computed at two thresholds (e.g. 1000 m for
  flash-flood QPE, 3000 m for cyclone structure) since a site can be
  excellent for one and poor for the other.
- **Largest contiguous blocked sector** — widest azimuthal run where
  `CBB > threshold`, treating azimuth as circular (a sector straddling due
  north is not split by an arbitrary array boundary). A small blockage
  fraction scattered in pinches is benign; the same fraction concentrated in
  one sector is a permanent hole.
- **Residual gap fraction** (network mode) — importance-weighted area seen
  by *no* radar in a joint configuration, at the higher threshold.

None of these are computed without an importance layer: an all-zero or fully
masked weight raster raises rather than silently returning a meaningless
`0/0` fraction (see `core/metrics.py:_weighted_fraction`).

## Nodata policy

`core/dem.py` never lets a nodata cell silently become `0` or `NaN` —  that
produces fictitious visibility that looks optimistic and plausible. Every
caller of `Dem.sample_lonlat` must choose a policy: `"raise"` (default, abort
on any nodata sample) or `"block"` (treat nodata as terrain that fully blocks
the beam, `+inf`). There is no silent default.

## Reproducibility manifest

Every algorithm run writes a manifest (`core/manifest.py`): timestamp,
plugin version, numpy/GDAL versions, every parameter that affects the result,
and a SHA-256 hash of every input raster/vector by content — not just by
path, since a path can be silently overwritten. It exists to close a single
objection — "how did you get this?" — for a report that justifies public
investment.
