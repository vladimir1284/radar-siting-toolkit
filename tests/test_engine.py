import numpy as np
import rasterio
from rasterio.transform import Affine

from core.dem import Dem
from core.engine import h_min_at_point, site_merit, sweep_site
from core.geometry import aeqd_proj4, transformer_from_aeqd


def _write_dem(path, terrain, res, lon0, lat0, nodata=-9999.0):
    n = terrain.shape[0]
    half_extent = n * res / 2.0
    transform = Affine(res, 0.0, -half_extent, 0.0, -res, half_extent)
    crs = aeqd_proj4(lon0, lat0)
    with rasterio.open(
        path,
        "w",
        driver="GTiff",
        height=n,
        width=n,
        count=1,
        dtype="float32",
        crs=crs,
        transform=transform,
        nodata=nodata,
    ) as dst:
        dst.write(terrain.astype("float32"), 1)
    return transform


def test_flat_terrain_is_fully_unblocked(tmp_path):
    """No obstacle anywhere: the beam center stays above flat ground at every
    azimuth and range, so CBB must be exactly 0 everywhere."""
    res = 100.0
    n = 200
    terrain = np.zeros((n, n), dtype="float32")
    path = tmp_path / "flat.tif"
    _write_dem(path, terrain, res, lon0=-76.8, lat0=18.0)

    with Dem(str(path)) as dem:
        out = sweep_site(
            dem,
            lon0=-76.8,
            lat0=18.0,
            h0=30.0,
            k=4 / 3,
            elevation_deg=0.5,
            beamwidth_deg=1.0,
            max_range_m=8_000.0,
            azimuth_step_deg=10.0,
        )

    assert out["azimuths_deg"].tolist() == list(np.arange(0.0, 360.0, 10.0))
    assert out["theta_rad"].shape == out["h_min_m"].shape == out["cbb"].shape
    assert out["theta_rad"].shape[0] == len(out["azimuths_deg"])
    assert np.all(out["cbb"] == 0.0)


def test_ridge_blocks_only_its_own_sector(tmp_path):
    """A ridge crossing azimuths 80-100 at ~4 km must drive CBB to 1 beyond
    it on those azimuths, and leave azimuth 0 (no ridge) fully clear."""
    res = 30.0
    n = 400
    half_extent = n * res / 2.0
    xs = -half_extent + res * (np.arange(n) + 0.5)
    ys = half_extent - res * (np.arange(n) + 0.5)
    X, Y = np.meshgrid(xs, ys)
    R = np.hypot(X, Y)
    AZ = np.degrees(np.arctan2(X, Y)) % 360.0

    terrain = np.zeros((n, n), dtype="float32")
    ridge_band = np.abs(R - 4000.0) < (res / 2.0)
    ridge_sector = (AZ > 80) & (AZ < 100)
    terrain = np.where(ridge_band & ridge_sector, 900.0, terrain)

    path = tmp_path / "ridge.tif"
    _write_dem(path, terrain, res, lon0=-76.8, lat0=18.0)

    with Dem(str(path)) as dem:
        out = sweep_site(
            dem,
            lon0=-76.8,
            lat0=18.0,
            h0=30.0,
            k=4 / 3,
            elevation_deg=0.5,
            beamwidth_deg=1.0,
            max_range_m=5_900.0,
            azimuth_step_deg=5.0,
        )

    az = out["azimuths_deg"]
    idx_through_ridge = np.argmin(np.abs(az - 90.0))
    idx_clear = np.argmin(np.abs(az - 0.0))

    assert out["cbb"][idx_through_ridge, -1] == 1.0
    assert out["cbb"][idx_clear, -1] == 0.0


def test_range_step_defaults_to_dem_native_resolution(tmp_path):
    res = 50.0
    n = 100
    terrain = np.zeros((n, n), dtype="float32")
    path = tmp_path / "flat.tif"
    _write_dem(path, terrain, res, lon0=-76.8, lat0=18.0)

    with Dem(str(path)) as dem:
        out = sweep_site(
            dem,
            lon0=-76.8,
            lat0=18.0,
            h0=10.0,
            k=4 / 3,
            elevation_deg=0.5,
            beamwidth_deg=1.0,
            max_range_m=2_000.0,
            azimuth_step_deg=90.0,
        )

    assert np.isclose(out["s_m"][0], res)
    assert np.allclose(np.diff(out["s_m"]), res)


def test_cbb_is_monotonic_non_decreasing_per_ray(tmp_path):
    res = 30.0
    n = 300
    rng = np.random.default_rng(3)
    terrain = rng.uniform(0, 400, size=(n, n)).astype("float32")
    path = tmp_path / "random.tif"
    _write_dem(path, terrain, res, lon0=-76.8, lat0=18.0)

    with Dem(str(path)) as dem:
        out = sweep_site(
            dem,
            lon0=-76.8,
            lat0=18.0,
            h0=50.0,
            k=4 / 3,
            elevation_deg=0.5,
            beamwidth_deg=1.0,
            max_range_m=4_000.0,
            azimuth_step_deg=15.0,
        )

    assert np.all(np.diff(out["cbb"], axis=1) >= -1e-9)
    assert np.all(np.diff(out["theta_rad"], axis=1) >= -1e-9)


def test_h_min_at_point_closer_than_range_step_returns_h0(tmp_path):
    res = 100.0
    n = 50
    terrain = np.zeros((n, n), dtype="float32")
    path = tmp_path / "flat.tif"
    lon0, lat0 = -76.8, 18.0
    _write_dem(path, terrain, res, lon0=lon0, lat0=lat0)
    to_lonlat = transformer_from_aeqd(lon0, lat0)

    with Dem(str(path)) as dem:
        query_lon, query_lat = to_lonlat.transform(10.0, 10.0)  # ~14 m away, < res
        h_min = h_min_at_point(
            dem, lon0, lat0, h0=30.0, k=4 / 3,
            query_lon=query_lon, query_lat=query_lat, range_step_m=res,
        )

    assert h_min == 30.0


def test_h_min_at_point_matches_only_its_own_bearing(tmp_path):
    """A ridge across azimuths 80-100 must raise h_min for a point behind it
    on that bearing (az=90), but not for an equally distant point on a clear
    bearing (az=0)."""
    res = 30.0
    n = 400
    half_extent = n * res / 2.0
    xs = -half_extent + res * (np.arange(n) + 0.5)
    ys = half_extent - res * (np.arange(n) + 0.5)
    X, Y = np.meshgrid(xs, ys)
    R = np.hypot(X, Y)
    AZ = np.degrees(np.arctan2(X, Y)) % 360.0

    terrain = np.zeros((n, n), dtype="float32")
    ridge_band = np.abs(R - 4000.0) < (res / 2.0)
    ridge_sector = (AZ > 80) & (AZ < 100)
    terrain = np.where(ridge_band & ridge_sector, 900.0, terrain)

    lon0, lat0 = -76.8, 18.0
    path = tmp_path / "ridge.tif"
    _write_dem(path, terrain, res, lon0=lon0, lat0=lat0)
    to_lonlat = transformer_from_aeqd(lon0, lat0)

    with Dem(str(path)) as dem:
        h0 = 30.0
        s_query = 5000.0
        lon_blocked, lat_blocked = to_lonlat.transform(
            s_query * np.sin(np.radians(90.0)), s_query * np.cos(np.radians(90.0))
        )
        lon_clear, lat_clear = to_lonlat.transform(
            s_query * np.sin(np.radians(0.0)), s_query * np.cos(np.radians(0.0))
        )
        h_min_blocked = h_min_at_point(
            dem, lon0, lat0, h0, k=4 / 3, query_lon=lon_blocked, query_lat=lat_blocked,
        )
        h_min_clear = h_min_at_point(
            dem, lon0, lat0, h0, k=4 / 3, query_lon=lon_clear, query_lat=lat_clear,
        )

    assert h_min_blocked > h_min_clear
    assert h_min_blocked > 100.0


def _ridge_dem_and_importance(tmp_path, lon0, lat0, importance_values):
    """Same ridge-across-80-100-degrees DEM as test_ridge_blocks_only_its_own_sector,
    plus an importance array on the identical grid (same resolution, extent,
    and AEQD CRS), so cell (i, j) means the same physical location in both."""
    res = 30.0
    n = 400
    half_extent = n * res / 2.0
    xs = -half_extent + res * (np.arange(n) + 0.5)
    ys = half_extent - res * (np.arange(n) + 0.5)
    X, Y = np.meshgrid(xs, ys)
    R = np.hypot(X, Y)
    AZ = np.degrees(np.arctan2(X, Y)) % 360.0

    terrain = np.zeros((n, n), dtype="float32")
    ridge_band = np.abs(R - 4000.0) < (res / 2.0)
    ridge_sector = (AZ > 80) & (AZ < 100)
    terrain = np.where(ridge_band & ridge_sector, 900.0, terrain)

    dem_path = tmp_path / "ridge.tif"
    transform = _write_dem(dem_path, terrain, res, lon0, lat0)
    importance = importance_values(AZ)
    return dem_path, transform, importance


def test_site_merit_drops_when_importance_is_behind_the_ridge(tmp_path):
    lon0, lat0 = -76.8, 18.0
    dem_path, transform, importance_blocked = _ridge_dem_and_importance(
        tmp_path, lon0, lat0, lambda az: np.where((az > 80) & (az < 100), 1.0, 0.0)
    )
    _, _, importance_clear = _ridge_dem_and_importance(
        tmp_path, lon0, lat0, lambda az: np.where((az > 80) & (az < 100), 0.0, 1.0)
    )
    crs = aeqd_proj4(lon0, lat0)

    with Dem(str(dem_path)) as dem:
        merit_blocked = site_merit(
            dem, importance_blocked, transform, crs, None,
            lon0, lat0, h0=30.0, k=4 / 3, elevation_deg=0.5, beamwidth_deg=1.0,
            max_range_m=5_900.0, azimuth_step_deg=5.0,
        )
        merit_clear = site_merit(
            dem, importance_clear, transform, crs, None,
            lon0, lat0, h0=30.0, k=4 / 3, elevation_deg=0.5, beamwidth_deg=1.0,
            max_range_m=5_900.0, azimuth_step_deg=5.0,
        )

    assert merit_blocked["merit"][1000.0] < merit_clear["merit"][1000.0]
    assert merit_clear["merit"][1000.0] > 0.9
    assert merit_blocked["largest_blocked_sector_deg"] == merit_clear["largest_blocked_sector_deg"]
