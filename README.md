# Radar Siting Toolkit

Open-source QGIS plugin for weather-radar coverage and siting analysis.

Given a territory, a set of constraints, and a definition of what area
matters, this toolkit answers: **where should one or more weather radars be
sited, and what will each configuration see?**

It is not a beam-blockage calculator. Partial beam blockage is one component
of the calculation, not the product. The product is reproducible geometric
evidence to support an investment recommendation to a government or a donor.

## Modes

- **Evaluate** — score a short list of candidate sites (or existing radars):
  coverage, blockage, and compared metrics per site and for the joint network.
- **Discover** — given an eligible-terrain mask, produce a merit raster over
  the whole territory to narrow the search space to a defensible short list.
- **Network** — joint coverage, redundancy, and residual gaps across a set of
  radars (own candidates plus existing domestic or foreign radars).

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

Implemented in plain numpy (no wradlib dependency at runtime — wradlib pulls
in scipy, matplotlib, xarray, xradar, and recommends conda install, which is
not installable from the QGIS Plugin Manager). wradlib is used only as an
offline validation reference, outside the plugin.

Projection is azimuthal equidistant (AEQD) centered on each radar, not UTM.

## Scope

Explicitly out of scope:

- Rain attenuation and frequency-band selection.
- Real anomalous propagation (ducting, measured refractivity gradients) —
  only the effective earth-radius factor is parameterized.
- Sea clutter, ground clutter, RF interference, spectrum compatibility.
- Quantitative precipitation estimation, nowcasting, radar data processing.
- Tower structure, radome, wind load, electrical resilience.
- Cost, lifecycle, staffing.

This toolkit covers roughly 15% of a full siting study.

## Elevation data

Recommended source: Copernicus GLO-30 (30 m, TanDEM-X derived). Because
GLO-30's interferometric phase center falls inside the canopy — neither a
faithful DSM nor a DTM — runs should be made with both GLO-30 (surface) and
FABDEM (bare earth) and the envelope between them reported as vegetation
uncertainty.

**FABDEM is CC BY-NC-SA 4.0 (non-commercial) and is never included in this
repository.** Resolve licensing in writing before using it in a paid
deliverable; without it, GLO-30 alone remains defensible with the limitation
declared.

## Status

Early scaffolding. Per the project plan, the analysis engine is developed and
validated in a notebook first (critical path for the originating field
report); the QGIS plugin wrapper (`processing/`, `core/`) follows afterward.
See `Plan_Tecnico_Radar_Siting_Toolkit.md` for the full technical plan and
work breakdown.

## License

GPL v2 or later (`LICENSE`) — the de facto requirement for QGIS plugins.
