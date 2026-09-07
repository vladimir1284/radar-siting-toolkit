"""Azimuthal-equidistant (AEQD) projection and ray sampling.

Plan_Tecnico_Radar_Siting_Toolkit.md, sections 3.1, 3.4, 3.5.

All geometry for a single radar's analysis happens in a local AEQD plane
centered on that radar: x/y in meters, distance from the radar is exact
Euclidean distance, and azimuth from the radar is exact bearing. This avoids
the anisotropic distortion UTM introduces at 250 km range and across zone
boundaries.
"""

import numpy as np
from pyproj import Geod, Transformer

EARTH_RADIUS_M = 6371000.0


def aeqd_proj4(lon0, lat0):
    """PROJ string for an AEQD projection centered on (lon0, lat0)."""
    return f"+proj=aeqd +lat_0={lat0} +lon_0={lon0} +datum=WGS84 +units=m +no_defs"


def transformer_to_aeqd(lon0, lat0):
    """Transformer from EPSG:4326 (lon, lat) to the local AEQD plane (x, y) in meters."""
    return Transformer.from_crs("EPSG:4326", aeqd_proj4(lon0, lat0), always_xy=True)


def transformer_from_aeqd(lon0, lat0):
    """Transformer from the local AEQD plane (x, y) in meters back to EPSG:4326 (lon, lat)."""
    return Transformer.from_crs(aeqd_proj4(lon0, lat0), "EPSG:4326", always_xy=True)


def ray_points(azimuth_deg, max_range_m, step_m):
    """Sample distances and local (x, y) offsets in meters along one ray from the origin.

    x is east, y is north, in the AEQD plane centered on the radar. Azimuth is
    degrees clockwise from north, matching standard radar convention.

    Returns
    -------
    s : (N,) ndarray of range along the ray [m], starting at step_m (not 0 —
        the radar's own cell is not part of the propagation path).
    x, y : (N,) ndarrays of local plane coordinates [m].
    """
    n = int(np.floor(max_range_m / step_m))
    s = (np.arange(1, n + 1)) * step_m
    az_rad = np.deg2rad(azimuth_deg)
    x = s * np.sin(az_rad)
    y = s * np.cos(az_rad)
    return s, x, y


def native_step(dem_transform, dem_crs, lat0=None):
    """Native ground sampling distance [m] of a DEM, from its affine transform.

    Section 3.5: sampling must use the DEM's own resolution, never a coarser
    fixed step, or intermediate ridges are skipped entirely.

    GLO-30 and FABDEM (config/jamaica.yaml) ship in EPSG:4326: dem_transform.a/.e
    are degrees, not meters. A geographic dem_crs converts its pixel size to
    meters via a geodesic step at lat0 (meters-per-degree-longitude shrinks
    with latitude, so this cannot be a fixed constant) -- treating a degrees
    transform as meters silently multiplies the ray point count by ~1e5.
    """
    dlon = abs(dem_transform.a)
    dlat = abs(dem_transform.e)
    if dem_crs is not None and dem_crs.is_geographic:
        if lat0 is None:
            raise ValueError(
                "native_step: dem_crs is geographic (degrees) -- lat0 is required "
                "to convert its pixel size to a ground distance in meters."
            )
        geod = Geod(ellps="WGS84")
        _, _, dx_m = geod.inv(0.0, lat0, dlon, lat0)
        _, _, dy_m = geod.inv(0.0, lat0, 0.0, lat0 + dlat)
        return min(dx_m, dy_m)
    return min(dlon, dlat)
