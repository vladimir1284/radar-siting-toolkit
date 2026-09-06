"""Importance-layer sampling along a radar sweep.

Plan_Tecnico_Radar_Siting_Toolkit.md, sections 4.1, 4.3, 6.2, 9 (Fase 3).

The primary figure of merit needs the importance-layer value at every point
of a candidate's own polar sweep (core.engine.sweep_site) -- the opposite
direction from alg_network (Fase 2), which resamples h_min onto the
importance layer's own Cartesian grid to join several radars over one
shared grid (core.engine.h_min_at_point). Sampling importance along a
single site's own rays is a plain raster read at each ray point, the same
O(N)-per-ray cost as the DEM sampling sweep_site already does.

Unlike core.dem.Dem, a nodata or out-of-extent importance cell here
contributes zero weight rather than raising or blocking (section 3.5's
fail-loud policy is specific to elevation, where silent zero produces
fictitious visibility; an unweighted importance cell is a modeling choice
-- no population there, or outside the area of interest -- not a
silent-bug risk).
"""

import numpy as np
import rasterio
from rasterio.warp import transform as warp_transform


def sample_importance_along_rays(band, transform, crs, nodata, lons, lats):
    """Importance value at each (lon, lat) point; 0.0 outside the raster
    extent or at a nodata cell.

    band : 2D ndarray, the importance raster's single band, already read.
    transform, crs, nodata : the same raster's affine transform, CRS, and
        declared nodata value (crs and nodata may be None).
    lons, lats : arrays of identical shape (e.g. sweep_site's per-ray
        (n_azimuths, n_ranges) lons_deg/lats_deg, in WGS84).
    """
    lons = np.asarray(lons, dtype=float)
    lats = np.asarray(lats, dtype=float)

    xs, ys = warp_transform("EPSG:4326", crs, lons.ravel().tolist(), lats.ravel().tolist())
    rows, cols = rasterio.transform.rowcol(transform, xs, ys)
    rows = np.asarray(rows)
    cols = np.asarray(cols)

    out_of_bounds = (rows < 0) | (rows >= band.shape[0]) | (cols < 0) | (cols >= band.shape[1])
    rows_clipped = np.clip(rows, 0, band.shape[0] - 1)
    cols_clipped = np.clip(cols, 0, band.shape[1] - 1)
    values = band[rows_clipped, cols_clipped].astype(float)

    invalid = out_of_bounds | np.isnan(values)
    if nodata is not None:
        invalid |= np.isclose(values, nodata)
    values = np.where(invalid, 0.0, values)

    return values.reshape(lons.shape)
