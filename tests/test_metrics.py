import numpy as np
import pytest

from core.metrics import (
    largest_contiguous_blocked_sector_deg,
    merit_scores,
    network_residual_gap_fraction,
    weighted_visible_fraction,
)


def test_weighted_visible_fraction_uniform_importance():
    h_min = np.array([500.0, 1500.0, 800.0, 4000.0])
    importance = np.ones_like(h_min)
    # 2 of 4 cells <= 1000 m
    assert weighted_visible_fraction(h_min, importance, 1000.0) == pytest.approx(0.5)


def test_weighted_visible_fraction_weights_by_importance():
    h_min = np.array([500.0, 5000.0])
    importance = np.array([9.0, 1.0])
    # visible cell carries 9x the weight of the blocked one
    assert weighted_visible_fraction(h_min, importance, 1000.0) == pytest.approx(0.9)


def test_weighted_visible_fraction_respects_mask():
    h_min = np.array([500.0, 5000.0, 500.0])
    importance = np.array([1.0, 1.0, 1.0])
    mask = np.array([True, True, False])
    assert weighted_visible_fraction(h_min, importance, 1000.0, mask=mask) == pytest.approx(0.5)


def test_zero_total_importance_raises():
    h_min = np.array([500.0, 500.0])
    importance = np.array([0.0, 0.0])
    with pytest.raises(ValueError):
        weighted_visible_fraction(h_min, importance, 1000.0)


def test_merit_scores_two_thresholds():
    h_min = np.array([500.0, 1500.0, 4000.0])
    importance = np.ones_like(h_min)
    scores = merit_scores(h_min, importance, thresholds=(1000.0, 3000.0))
    assert scores[1000.0] == pytest.approx(1 / 3)
    assert scores[3000.0] == pytest.approx(2 / 3)


def test_largest_contiguous_blocked_sector_simple_run():
    cbb = np.zeros(360)
    cbb[30:100] = 0.9  # blocked from 30 to 99 inclusive -> 70 deg
    assert largest_contiguous_blocked_sector_deg(cbb, azimuth_step_deg=1.0) == pytest.approx(70.0)


def test_largest_contiguous_blocked_sector_wraps_across_zero():
    cbb = np.zeros(360)
    cbb[350:] = 0.9
    cbb[:10] = 0.9  # 10 before + 10 after wrap = 20 deg contiguous through 0
    assert largest_contiguous_blocked_sector_deg(cbb, azimuth_step_deg=1.0) == pytest.approx(20.0)


def test_largest_contiguous_blocked_sector_all_clear():
    cbb = np.full(360, 0.1)
    assert largest_contiguous_blocked_sector_deg(cbb, azimuth_step_deg=1.0) == 0.0


def test_largest_contiguous_blocked_sector_fully_blocked():
    cbb = np.full(360, 0.9)
    assert largest_contiguous_blocked_sector_deg(cbb, azimuth_step_deg=1.0) == pytest.approx(360.0)


def test_network_residual_gap_needs_every_radar_to_miss():
    # radar 0 sees cell 1 well, radar 1 sees cell 0 well -> no residual gap
    h_min_stack = np.array([
        [500.0, 4000.0],
        [4000.0, 500.0],
    ])
    importance = np.array([1.0, 1.0])
    assert network_residual_gap_fraction(h_min_stack, importance, threshold_m=3000.0) == 0.0


def test_network_residual_gap_when_all_radars_miss_a_cell():
    h_min_stack = np.array([
        [500.0, 4000.0],
        [500.0, 5000.0],
    ])
    importance = np.array([1.0, 1.0])
    assert network_residual_gap_fraction(h_min_stack, importance, threshold_m=3000.0) == pytest.approx(0.5)
