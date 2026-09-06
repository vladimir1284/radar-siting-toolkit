import numpy as np
import pytest
import rasterio
from rasterio.transform import from_origin, xy as transform_xy

from core.dem import Dem, DemNodataError

NODATA_VALUE = -9999.0


def _write_synthetic_dem(path, values, nodata=NODATA_VALUE, crs="EPSG:4326"):
    """4x4 grid, 1 degree cells, north-up, origin at (lon=0, lat=4)."""
    transform = from_origin(0, 4, 1, 1)
    with rasterio.open(
        path,
        "w",
        driver="GTiff",
        height=values.shape[0],
        width=values.shape[1],
        count=1,
        dtype=values.dtype,
        crs=crs,
        transform=transform,
        nodata=nodata,
    ) as dst:
        dst.write(values, 1)
    return transform


def _cell_lonlat(transform, row, col):
    """Center of pixel (row, col), to sample exactly that cell back."""
    lon, lat = transform_xy(transform, row, col)
    return lon, lat


def test_sample_lonlat_reads_correct_elevation(tmp_path):
    values = np.array(
        [
            [10.0, 20.0, 30.0, 40.0],
            [50.0, 60.0, 70.0, 80.0],
            [90.0, 100.0, 110.0, 120.0],
            [130.0, 140.0, 150.0, 160.0],
        ],
        dtype="float32",
    )
    path = tmp_path / "dem.tif"
    transform = _write_synthetic_dem(path, values)

    with Dem(str(path)) as dem:
        lon, lat = _cell_lonlat(transform, 1, 2)
        z = dem.sample_lonlat([lon], [lat])
        assert z[0] == pytest.approx(70.0)


def test_declared_nodata_value_raises_by_default(tmp_path):
    values = np.full((4, 4), 100.0, dtype="float32")
    values[2, 1] = NODATA_VALUE
    path = tmp_path / "dem.tif"
    transform = _write_synthetic_dem(path, values)

    with Dem(str(path)) as dem:
        lon, lat = _cell_lonlat(transform, 2, 1)
        with pytest.raises(DemNodataError):
            dem.sample_lonlat([lon], [lat])


def test_nan_pixel_is_nodata_even_without_declared_nodata(tmp_path):
    values = np.full((4, 4), 100.0, dtype="float32")
    values[0, 0] = np.nan
    path = tmp_path / "dem.tif"
    transform = _write_synthetic_dem(path, values, nodata=None)

    with Dem(str(path)) as dem:
        assert dem.nodata is None
        lon, lat = _cell_lonlat(transform, 0, 0)
        with pytest.raises(DemNodataError):
            dem.sample_lonlat([lon], [lat])


def test_nodata_never_silently_becomes_zero_or_nan(tmp_path):
    """Section 3.5: the one bug that must be impossible, checked directly."""
    values = np.full((4, 4), 100.0, dtype="float32")
    values[3, 3] = NODATA_VALUE
    path = tmp_path / "dem.tif"
    transform = _write_synthetic_dem(path, values)

    with Dem(str(path)) as dem:
        lon, lat = _cell_lonlat(transform, 3, 3)
        z = dem.sample_lonlat([lon], [lat], nodata_policy="block")
        assert np.isinf(z[0]) and z[0] > 0
        assert z[0] != 0.0
        assert not np.isnan(z[0])


def test_block_policy_leaves_valid_cells_untouched(tmp_path):
    values = np.full((4, 4), 100.0, dtype="float32")
    values[3, 3] = NODATA_VALUE
    path = tmp_path / "dem.tif"
    transform = _write_synthetic_dem(path, values)

    with Dem(str(path)) as dem:
        lon0, lat0 = _cell_lonlat(transform, 0, 0)
        lon1, lat1 = _cell_lonlat(transform, 3, 3)
        z = dem.sample_lonlat([lon0, lon1], [lat0, lat1], nodata_policy="block")
        assert z[0] == pytest.approx(100.0)
        assert np.isinf(z[1])


def test_unknown_nodata_policy_raises(tmp_path):
    """Only exercised on a ray that actually hits nodata: with none, the value
    passes through untouched regardless of policy string, per section 3.5 --
    "block"/"raise" only ever apply to the nodata cells themselves."""
    values = np.full((4, 4), 100.0, dtype="float32")
    values[1, 1] = NODATA_VALUE
    path = tmp_path / "dem.tif"
    transform = _write_synthetic_dem(path, values)

    with Dem(str(path)) as dem:
        lon, lat = _cell_lonlat(transform, 1, 1)
        with pytest.raises(ValueError):
            dem.sample_lonlat([lon], [lat], nodata_policy="bogus")


def test_out_of_bounds_point_raises(tmp_path):
    values = np.full((4, 4), 100.0, dtype="float32")
    path = tmp_path / "dem.tif"
    _write_synthetic_dem(path, values)

    with Dem(str(path)) as dem:
        with pytest.raises(DemNodataError):
            dem.sample_lonlat([500.0], [500.0])


def test_multi_band_raster_rejected(tmp_path):
    path = tmp_path / "dem.tif"
    transform = from_origin(0, 4, 1, 1)
    with rasterio.open(
        path,
        "w",
        driver="GTiff",
        height=4,
        width=4,
        count=2,
        dtype="float32",
        crs="EPSG:4326",
        transform=transform,
    ) as dst:
        dst.write(np.zeros((4, 4), dtype="float32"), 1)
        dst.write(np.zeros((4, 4), dtype="float32"), 2)

    with pytest.raises(ValueError):
        Dem(str(path))
