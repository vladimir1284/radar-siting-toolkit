import numpy as np
from pyproj import Geod

from core.geometry import ray_points, transformer_to_aeqd, transformer_from_aeqd


def test_ray_points_length_and_range_convention():
    s, x, y = ray_points(azimuth_deg=0, max_range_m=1000, step_m=250)
    assert len(s) == 4
    assert s[0] == 250  # radar's own cell is excluded
    assert s[-1] == 1000


def test_ray_points_azimuth_directions():
    # North: x (east) ~ 0, y (north) > 0
    _, x, y = ray_points(0, 1000, 500)
    assert np.allclose(x, 0, atol=1e-9)
    assert np.all(y > 0)
    # East: x > 0, y ~ 0
    _, x, y = ray_points(90, 1000, 500)
    assert np.allclose(y, 0, atol=1e-9)
    assert np.all(x > 0)


def test_aeqd_round_trip():
    lon0, lat0 = -77.3, 18.0  # Kingston, Jamaica, approx.
    to_local = transformer_to_aeqd(lon0, lat0)
    to_lonlat = transformer_from_aeqd(lon0, lat0)

    x, y = to_local.transform(-76.9, 18.3)
    lon, lat = to_lonlat.transform(x, y)
    assert np.isclose(lon, -76.9, atol=1e-6)
    assert np.isclose(lat, 18.3, atol=1e-6)


def test_aeqd_distance_matches_range():
    """Distance from the origin in the AEQD plane must equal geodesic range."""
    lon0, lat0 = -77.3, 18.0
    to_local = transformer_to_aeqd(lon0, lat0)
    x, y = to_local.transform(lon0, lat0 + 0.1)  # ~11.1 km due north

    geod = Geod(ellps="WGS84")
    _, _, geodesic_dist = geod.inv(lon0, lat0, lon0, lat0 + 0.1)

    r = np.hypot(x, y)
    assert np.isclose(r, geodesic_dist, rtol=1e-3)
