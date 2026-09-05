"""Accumulated terrain-horizon angle, and minimum visible height.

Plan_Tecnico_Radar_Siting_Toolkit.md, section 3.1.

One O(N) pass per ray, replacing the discrete elevation-angle loop of the
original plan: for each range s, compute the apparent terrain angle above
the radar (corrected for earth curvature via the effective earth-radius
factor k), then take the running maximum. Everything else (minimum clear
elevation, minimum visible height) is derived directly from that one array.
"""

import numpy as np

from .geometry import EARTH_RADIUS_M


def apparent_terrain_height(z, s, h0, k, earth_radius=EARTH_RADIUS_M):
    """Delta z(s): terrain height relative to the radar, net of earth curvature drop.

    z, s : (N,) ndarrays, terrain elevation [m] and range [m] along the ray.
    h0 : radar antenna height [m] = ground elevation at the radar + tower height.
    k : effective earth-radius factor (1.0, 4/3, or a superrefractive value).
    """
    return z - h0 - (s**2) / (2.0 * k * earth_radius)


def horizon_angle(z, s, h0, k, earth_radius=EARTH_RADIUS_M):
    """Theta(s): cumulative (running-max) horizon angle [rad] along the ray.

    This is the minimum clear elevation angle at each range: below it, the
    line of sight to that range is blocked by nearer terrain.
    """
    dz = apparent_terrain_height(z, s, h0, k, earth_radius)
    theta_terrain = np.arctan2(dz, s)
    return np.maximum.accumulate(theta_terrain)


def min_visible_height(z, s, h0, k, theta=None, earth_radius=EARTH_RADIUS_M):
    """H_min(s): minimum height [m] above local terrain visible from the radar at range s.

    If `theta` (from `horizon_angle`) is not supplied, it is computed here.
    """
    if theta is None:
        theta = horizon_angle(z, s, h0, k, earth_radius)
    return h0 + s * np.tan(theta) + (s**2) / (2.0 * k * earth_radius) - z
