# QGIS Algorithms

All four algorithms are registered under **Processing Toolbox → Radar Siting
Toolkit** (`processing/provider.py`). Each wraps a `core/` engine call for
QGIS I/O and writes an execution [manifest](core-reference.md#coremanifest)
alongside its outputs.

## Evaluate

`processing/alg_evaluate.py` — Fase 1 of the plan.

Runs the horizon-angle and Bech engine (`core.engine.sweep_site`) over a
full 360° scan from each candidate site (or existing radar), against a
single DEM.

**Input:** point layer of candidate sites, one DEM.

**Key fields expected on the input layer:** `tower_height_m` (required),
optional `ground_elevation_m` (sampled from the DEM if absent), optional
`name` (falls back to `site_<fid>`).

**Output:**

- Site summary layer: input points + `ground_elevation_m`, `h0_m`,
  `largest_blocked_sector_deg`, `worst_az_min_clear_elev_deg`,
  `worst_az_cbb_at_max_range`, `azimuth_csv` (path).
- One per-azimuth CSV per site: `azimuth_deg`,
  `min_clear_elevation_deg_at_max_range`, `h_min_m_at_max_range`,
  `cbb_at_max_range`.
- `manifest.json` for the whole run.

**Not included:** importance-weighted merit scores (needs an
area-of-interest raster grid, not just the per-site ray sweep) — see
[Robustness matrix](#robustness-matrix).

## Discover

`processing/alg_discover.py` — Fase 2 of the plan.

Runs a full 360° engine sweep from *every* eligible cell of a terrain mask,
with one uniform assumed tower height, and writes a merit raster on the
mask's own grid. **Does not decide a site** — narrows a whole territory to a
defensible short list for Evaluate.

**Input:** DEM, eligible-terrain mask (nonzero = eligible).

**Output:** 3-band raster — area-weighted visible fraction at two height
thresholds, and largest contiguous blocked sector, per eligible cell — plus
a manifest.

Merit here is **area-weighted** (`s·ds·dθ`, since equal azimuth/range steps
cover more physical area at longer range), not importance-weighted: this
mode's grid is a search space (every eligible cell probed as a hypothetical
site), not the area of interest an importance layer describes. Coarser
default azimuth step (5°) than Evaluate (1°) trades accuracy for the larger
cell count.

## Network

`processing/alg_network.py` — Fase 2 of the plan.

Joint coverage across a set of radars (own candidates plus existing radars,
domestic or foreign), over the importance layer's own grid. For every
radar/cell pair within range, computes the **exact** minimum visible height
(`core.engine.h_min_at_point` — one ray along the exact bearing, not an
interpolation between two of `sweep_site`'s fixed azimuths, since a cell
essentially never lands on one), then takes the per-cell minimum across
radars.

**Input:** point layer of radars, DEM, importance layer.

**Output:**

- Raster of best (minimum) `h_min` across all radars.
- `..._summary.json`: per-radar label/position/`h0`, joint visible fraction
  at both thresholds, residual gap fraction at the higher threshold.
- Manifest.

Unlike Discover, this uses the importance layer's real weights — the output
is the importance-weighted figure of merit for a *given* network
configuration, not a search-space proxy.

## Robustness matrix

`processing/alg_robustness.py` — Fase 3 of the plan.

For every candidate site, sweeps the importance-weighted figure of merit
(`core.engine.site_merit`) across every combination of DEM source,
effective earth-radius factor `k`, and tower height, then ranks candidates
within each combination.

**Input:** candidate sites, one or more DEM sources to sweep (e.g. GLO-30
and FABDEM), importance layer, comma-separated `k` values (default
`1.0,1.3333,2.0`), comma-separated tower heights (default `10,15,20,25,30`).

Ranking is driven by the **higher** height threshold (threshold 2, e.g. 3 km
— typically the more demanding "can it see cyclone structure" number); both
thresholds are reported per row so a reviewer can re-rank by the other one
by hand.

**Output:**

- `robustness_matrix.csv` — one row per (site, DEM source, k, tower height):
  merit at both thresholds, largest blocked sector, rank within that combo.
- `marginal_performance.csv` — merit gained per additional meter of tower
  height, for each (site, DEM source, k) line — the budget argument.
- `summary.json` — whether the ranking is identical across every
  combination (`ranking_stable_across_all_combinations`), the single top
  site if there is one, and each site's rank range.
- `figures/` (optional, requires matplotlib, on by default): merit-vs-tower
  plot per site, plus a rank-stability bar chart.
- Manifest.

If the ranking flips across combinations, the team finds out before a
reviewer does — this is the mechanism section 6.2 of the plan calls
"shielding the recommendation." `matplotlib` is a soft dependency, lazily
imported: CSV/JSON output is unaffected if it's missing.
