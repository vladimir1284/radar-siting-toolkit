# Radar Siting Toolkit

Open-source QGIS plugin for weather-radar coverage and siting analysis.

Given a territory, a set of constraints, and a definition of what area
matters, this toolkit answers: **where should one or more weather radars be
sited, and what will each configuration see?**

It is not a beam-blockage calculator. Partial beam blockage is one component
of the calculation, not the product. The product is reproducible geometric
evidence to support an investment recommendation to a government or a donor —
every run writes a [manifest](core-reference.md#coremanifest) proving which
inputs, parameters, and library versions produced it.

This toolkit covers roughly **15%** of a full siting study. See
[Scope & Status](status.md) for what is deliberately excluded.

## Modes

| Mode | Input | Output | Use when |
|---|---|---|---|
| [Evaluate](algorithms.md#evaluate) | Point layer of candidate sites | Per-site coverage/blockage summary + manifest | You already have a short list |
| [Discover](algorithms.md#discover) | Eligible-terrain mask | Merit raster over the whole search space | You need to narrow a territory to a short list |
| [Network](algorithms.md#network) | Set of radars + importance layer | Joint coverage, redundancy, residual gaps | You're combining candidates with existing radars |
| [Robustness matrix](algorithms.md#robustness-matrix) | Candidate sites + DEM sources + parameter ranges | Ranking stability across k / tower height / DEM | You need to defend a recommendation against assumption changes |

## Core method

For each azimuth, the DEM is sampled along the ray at native resolution and
the maximum horizon angle seen up to each range is accumulated, corrected for
earth curvature via an effective earth-radius factor `k`. From this, minimum
clear elevation and minimum visible height are derived directly, in a single
`O(N)` pass per ray — no discrete elevation-angle loop.

For short-list candidates, partial beam blockage uses the circular-segment
geometry of Bech et al. (2003), always published as the range-cumulative
maximum (CBB), never as raw per-bin PBB — publishing per-bin PBB produces
wrong maps across the whole sector behind an obstacle.

See [Methodology](methodology.md) for the full derivation and
[Core API Reference](core-reference.md) for the implementation.

Implemented in plain numpy — no `wradlib` dependency at runtime (`wradlib`
pulls in scipy, matplotlib, xarray, xradar, and recommends a conda install,
which is not installable from the QGIS Plugin Manager). `wradlib` is used
only as an offline validation reference, outside the plugin.

Projection is azimuthal equidistant (AEQD) centered on each radar, not UTM —
see [Methodology → Geometry](methodology.md#geometry-aeqd-projection).

## Elevation data

Recommended source: Copernicus GLO-30 (30 m, TanDEM-X derived). Because
GLO-30's interferometric phase center falls inside the canopy — neither a
faithful DSM nor a DTM — runs should be made with both GLO-30 (surface) and
FABDEM (bare earth) and the envelope between them reported as vegetation
uncertainty. See [Elevation Data & Configuration](data.md).

## License

GPL v2 or later (`LICENSE`) — the de facto requirement for QGIS plugins.
