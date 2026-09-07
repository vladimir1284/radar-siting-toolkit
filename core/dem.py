"""DEM access with explicit nodata handling.

Plan_Tecnico_Radar_Siting_Toolkit.md, sections 3.5, 5.5.

Section 3.5 names this the most likely and hardest-to-detect bug in the
engine: a nodata cell that silently becomes 0 or NaN produces fictitious
visibility that looks optimistic and plausible. There is no implicit
default here — every caller must choose a policy, and "raise" is what you
get if you don't.
"""

import numpy as np
import rasterio
from rasterio.warp import transform as warp_transform


class DemNodataError(ValueError):
    """Raised when a sampled ray crosses nodata and the policy is 'raise'."""


class Dem:
    """A single-band elevation raster opened for point sampling.

    Parameters
    ----------
    path : str
        Path to a GDAL-readable raster (e.g. GLO-30 or FABDEM GeoTIFF).
    vertical_datum : str
        Free-text label recorded for the manifest (section 5.5), e.g.
        "EGM2008 geoid" for GLO-30. Not used in computation.
    """

    def __init__(self, path, vertical_datum=None):
        self.path = path
        self.vertical_datum = vertical_datum
        self.dataset = rasterio.open(path)
        if self.dataset.count != 1:
            raise ValueError(f"{path}: expected a single-band elevation raster, got {self.dataset.count} bands")
        self.nodata = self.dataset.nodata
        # Read once: sample_lonlat is called once per azimuth (hundreds of times
        # per sweep_site), and re-reading the full band from disk each call
        # made a single site's sweep re-read a ~100 MB array hundreds of times.
        self._band = self.dataset.read(1)

    @property
    def transform(self):
        return self.dataset.transform

    @property
    def crs(self):
        return self.dataset.crs

    def sample_lonlat(self, lons, lats, nodata_policy="raise"):
        """Sample elevation at WGS84 (lon, lat) points, reprojecting to the DEM's CRS.

        Parameters
        ----------
        nodata_policy : {"raise", "block"}
            "raise" (default): raise DemNodataError if any sampled cell is
            nodata. "block": return +inf for nodata cells, i.e. treat them as
            terrain that fully blocks the beam, per section 3.5. Never
            silently returns 0 or NaN for a nodata cell.
        """
        lons = np.asarray(lons, dtype=float)
        lats = np.asarray(lats, dtype=float)
        xs, ys = warp_transform("EPSG:4326", self.crs, lons.tolist(), lats.tolist())
        return self._sample_xy(np.asarray(xs), np.asarray(ys), nodata_policy)

    def _sample_xy(self, xs, ys, nodata_policy):
        rows, cols = rasterio.transform.rowcol(self.transform, xs, ys)
        rows = np.asarray(rows)
        cols = np.asarray(cols)

        band = self._band
        out_of_bounds = (
            (rows < 0) | (rows >= band.shape[0]) | (cols < 0) | (cols >= band.shape[1])
        )
        if out_of_bounds.any():
            raise DemNodataError(
                f"{out_of_bounds.sum()} sampled point(s) fall outside the DEM extent "
                f"({self.path}); widen the input raster before analysis."
            )

        z = band[rows, cols].astype(float)
        is_nodata = np.zeros_like(z, dtype=bool)
        if self.nodata is not None:
            is_nodata |= np.isclose(z, self.nodata)
        is_nodata |= np.isnan(z)

        if is_nodata.any():
            if nodata_policy == "raise":
                raise DemNodataError(
                    f"{int(is_nodata.sum())} sampled point(s) hit nodata in {self.path}. "
                    "Fill or clip the DEM, or pass nodata_policy='block' to treat these "
                    "cells as fully blocking terrain."
                )
            elif nodata_policy == "block":
                z = np.where(is_nodata, np.inf, z)
            else:
                raise ValueError(f"unknown nodata_policy: {nodata_policy!r}")

        return z

    def close(self):
        self.dataset.close()

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()


def raster_cell_centers_lonlat(transform, crs, height, width):
    """WGS84 (lon, lat) of every cell center of a raster grid.

    Used by alg_discover / alg_network (Fase 2) to turn a mask or importance
    raster's own grid into query points for the radar engine, which works in
    WGS84 lon/lat throughout (section 3.4). Returns two (height, width)
    arrays, lons and lats.
    """
    rows, cols = np.indices((height, width))
    xs, ys = rasterio.transform.xy(transform, rows.ravel(), cols.ravel(), offset="center")
    lons, lats = warp_transform(crs, "EPSG:4326", xs, ys)
    lons = np.asarray(lons).reshape(height, width)
    lats = np.asarray(lats).reshape(height, width)
    return lons, lats
