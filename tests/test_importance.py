import numpy as np
import pytest
from rasterio.transform import from_origin, xy as transform_xy

from core.importance import sample_importance_along_rays


def _cell_lonlat(transform, row, col):
    lon, lat = transform_xy(transform, row, col)
    return lon, lat


def test_samples_known_cell_values():
    values = np.array(
        [
            [10.0, 20.0, 30.0, 40.0],
            [50.0, 60.0, 70.0, 80.0],
            [90.0, 100.0, 110.0, 120.0],
            [130.0, 140.0, 150.0, 160.0],
        ]
    )
    transform = from_origin(0, 4, 1, 1)
    crs = "EPSG:4326"

    lon, lat = _cell_lonlat(transform, 1, 2)
    out = sample_importance_along_rays(values, transform, crs, None, [lon], [lat])
    assert out[0] == pytest.approx(70.0)


def test_nodata_cell_contributes_zero_not_raise():
    values = np.array(
        [
            [10.0, -9999.0],
            [30.0, 40.0],
        ]
    )
    transform = from_origin(0, 2, 1, 1)
    crs = "EPSG:4326"

    lon, lat = _cell_lonlat(transform, 0, 1)
    out = sample_importance_along_rays(values, transform, crs, -9999.0, [lon], [lat])
    assert out[0] == 0.0


def test_out_of_bounds_point_contributes_zero_not_raise():
    values = np.full((4, 4), 100.0)
    transform = from_origin(0, 4, 1, 1)
    crs = "EPSG:4326"

    out = sample_importance_along_rays(values, transform, crs, None, [500.0], [500.0])
    assert out[0] == 0.0


def test_nan_cell_contributes_zero():
    values = np.array([[np.nan, 10.0], [20.0, 30.0]])
    transform = from_origin(0, 2, 1, 1)
    crs = "EPSG:4326"

    lon, lat = _cell_lonlat(transform, 0, 0)
    out = sample_importance_along_rays(values, transform, crs, None, [lon], [lat])
    assert out[0] == 0.0


def test_preserves_input_shape():
    values = np.full((10, 10), 5.0)
    transform = from_origin(0, 10, 1, 1)
    crs = "EPSG:4326"

    lons = np.array([[0.5, 1.5], [2.5, 3.5]])
    lats = np.array([[9.5, 8.5], [7.5, 6.5]])
    out = sample_importance_along_rays(values, transform, crs, None, lons, lats)
    assert out.shape == (2, 2)
    assert np.all(out == 5.0)
