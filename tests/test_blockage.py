import numpy as np
import pytest

from core.blockage import beam_block_frac, cum_beam_block_frac


def test_boundary_values():
    a = 100.0
    assert beam_block_frac(np.array([0.0]), np.array([a]))[0] == pytest.approx(0.5)
    assert beam_block_frac(np.array([a]), np.array([a]))[0] == pytest.approx(1.0)
    assert beam_block_frac(np.array([-a]), np.array([a]))[0] == pytest.approx(0.0)


def test_monotonic_in_y():
    a = np.full(5, 50.0)
    y = np.array([-60.0, -25.0, 0.0, 25.0, 60.0])
    pbb = beam_block_frac(y, a)
    assert np.all(np.diff(pbb) >= 0)


def test_cumulative_is_running_max():
    pbb = np.array([[0.1, 0.9, 0.2, 0.05, 0.3]])
    cbb = cum_beam_block_frac(pbb, axis=1)
    assert np.array_equal(cbb[0], [0.1, 0.9, 0.9, 0.9, 0.9])


def test_cumulative_never_decreases_along_range():
    rng = np.random.default_rng(0)
    pbb = rng.uniform(0, 1, size=(20, 200))
    cbb = cum_beam_block_frac(pbb, axis=1)
    assert np.all(np.diff(cbb, axis=1) >= -1e-12)


def test_matches_wradlib():
    wradlib = pytest.importorskip("wradlib")
    rng = np.random.default_rng(1)
    th = rng.uniform(-200, 200, size=1000)
    bh = np.zeros_like(th)
    a = np.full_like(th, 80.0)

    ours = beam_block_frac(th - bh, a)
    reference = wradlib.qual.beam_block_frac(th, bh, a)
    np.testing.assert_allclose(ours, reference, atol=1e-9)
