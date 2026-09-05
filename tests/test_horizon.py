import numpy as np

from core.horizon import horizon_angle, min_visible_height


def test_flat_terrain_small_angle_near_range():
    """At short range over flat ground, theta ~ -h0/s (small-angle, curvature negligible)."""
    s = np.array([100.0])
    z = np.zeros_like(s)
    h0 = 30.0
    theta = horizon_angle(z, s, h0, k=4 / 3)
    assert np.isclose(theta[0], np.arctan(-h0 / s[0]), rtol=1e-3)


def test_running_max_never_decreases():
    rng = np.random.default_rng(2)
    s = np.linspace(30, 100_000, 2000)
    z = rng.uniform(0, 500, size=s.shape)
    theta = horizon_angle(z, s, h0=25.0, k=4 / 3)
    assert np.all(np.diff(theta) >= -1e-12)


def test_ridge_creates_a_shadow():
    """A single ridge must raise the horizon angle for every range beyond it."""
    s = np.linspace(100, 20_000, 400)
    z = np.zeros_like(s)
    ridge_idx = 50
    z[ridge_idx] = 800.0  # a sharp 800 m ridge
    h0 = 20.0
    theta = horizon_angle(z, s, h0, k=4 / 3)
    # theta at the ridge and everywhere beyond must equal theta at the ridge
    assert np.allclose(theta[ridge_idx:], theta[ridge_idx])
    # and it must exceed theta just before the ridge
    assert theta[ridge_idx] > theta[ridge_idx - 1]


def test_min_visible_height_is_zero_at_the_horizon_and_nonnegative_elsewhere():
    """H_min(s) is height above *local terrain*, not absolute height: by
    construction it is exactly 0 at whichever range currently sets the
    running-max horizon (the line of sight grazes that point's terrain),
    and non-negative everywhere (never claims sub-ground visibility)."""
    s = np.linspace(100, 20_000, 400)
    z = np.zeros_like(s)
    z[50] = 800.0
    h0 = 20.0
    theta = horizon_angle(z, s, h0, k=4 / 3)
    hmin = min_visible_height(z, s, h0, k=4 / 3, theta=theta)
    # The ridge itself sets the horizon: line of sight exactly grazes it.
    assert np.isclose(hmin[50], 0.0, atol=1e-6)
    assert np.all(hmin >= -1e-6)
