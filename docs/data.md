# Elevation Data & Configuration

## Recommended sources

**Copernicus GLO-30** (COP30, 30 m, TanDEM-X derived) is the recommended
default. Its interferometric phase center falls *inside the canopy* — it is
neither a faithful DSM (surface, including canopy) nor a faithful DTM (bare
earth). Runs should therefore be made with both:

- **GLO-30** — surface-ish, via the OpenTopography `globaldem` API.
- **FABDEM** — bare-earth correction of GLO-30/GLO-90.

and the envelope between them reported as vegetation uncertainty. The
[Robustness matrix](algorithms.md#robustness-matrix) algorithm sweeps both
sources automatically for exactly this reason.

!!! warning "FABDEM licensing"
    FABDEM is **CC BY-NC-SA 4.0 (non-commercial)** and must never be
    committed to this repository. Resolve licensing in writing before using
    it in a paid deliverable; without it, GLO-30 alone remains defensible
    with the limitation declared. `config/jamaica.yaml`'s FABDEM entry is
    marked dev/validation-only for this reason.

## DEM vertical datum

Both GLO-30 and FABDEM report elevation relative to the **EGM2008 geoid**.
The `DEM_VERTICAL_DATUM` parameter on every algorithm is recorded in the
manifest for traceability only — it is not used in computation, since the
engine works entirely in relative heights (`H_min`, beam height) referenced
to the radar's own antenna height `h0`.

## Case configuration format

`config/*.yaml` describes one case's radars, DEM sources, importance layer,
and scan parameters. See `config/jamaica.yaml`:

```yaml
radars:
  - name: "Cooper's Hill Doppler Weather Radar, Saint Andrew Parish"
    lat: 18.07330822439374
    lon: -76.8527864446781
    ground_elevation_m: 765.6  # FABDEM bare-earth sample at this point
    tower_height_m: 30  # measured antenna height AGL, existing structure
dem:
  sources:
    - path: examples/jamaica/glo30_jamaica.tif
      label: GLO-30
      vertical_datum: EGM2008 geoid
    - path: examples/jamaica/fabdem_jamaica.tif
      label: FABDEM
      vertical_datum: EGM2008 geoid
importance_layer: null
scan:
  earth_radius_factor_k: [1.0, 1.3333, 2.0]  # normal / 4/3 / superrefractive
  tower_height_m: [10, 15, 20, 25, 30]
```

This file is a plan/example artifact for the Jamaica case, not yet consumed
by the QGIS algorithms directly (each algorithm takes its parameters through
the Processing dialog) — see `notebooks/fase0_motor_horizonte_bech.ipynb`
for where it is read today.

## Example data

`examples/jamaica/` ships cropped GLO-30 and FABDEM rasters over Jamaica
(south=17.65, north=18.55, west=-78.4, east=-76.2) plus a sample
`manifest.json`, for running the notebook and the plugin against a real
case without downloading data first. This extent does not cover Cuba/Haiti
terrain within a 250 km radar range — traded off to keep the repository
checkout small.

## Effective earth-radius factor `k`

Standard values used across the toolkit's defaults and `config/jamaica.yaml`:

| `k` | Meaning |
|---|---|
| `1.0` | No atmospheric refraction correction |
| `4/3` (`1.3333`) | Standard tropospheric refraction (the conventional default) |
| `2.0` | Superrefractive conditions |

`k` is swept, not fixed, in the Robustness matrix — the tropics can depart
from the standard 4/3 assumption, and the plan treats this as a modeling
uncertainty to be shown, not resolved by picking one value.
