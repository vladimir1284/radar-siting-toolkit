"""Partial beam blockage (Bech et al. 2003), reimplemented in plain numpy.

Plan_Tecnico_Radar_Siting_Toolkit.md, sections 3.2, 3.3.

This is the ~50-line reimplementation the plan calls for, dropping the
wradlib runtime dependency (scipy + matplotlib + xarray + xradar, conda-only
install path — not installable from the QGIS Plugin Manager). wradlib is
used only offline, in this notebook, to validate this module (section 7,
Nivel 2).

The cumulative form (`cum_beam_block_frac`) is mandatory for any published
map: publishing per-bin PBB instead produces wrong results across the whole
sector behind an obstacle, because a beam once blocked stays blocked at
every greater range regardless of what the terrain does next.
"""

import numpy as np

from .geometry import EARTH_RADIUS_M


def beam_height(s, h0, elevation_rad, k, earth_radius=EARTH_RADIUS_M):
    """Height [m] of the beam center at range s, corrected for earth curvature."""
    return h0 + s * np.tan(elevation_rad) + (s**2) / (2.0 * k * earth_radius)


def half_power_radius(s, beamwidth_rad):
    """a(s): half-power (-3 dB) beam radius [m] at range s."""
    return s * np.tan(beamwidth_rad / 2.0)


def beam_block_frac(y, a):
    """PBB: partial beam blockage fraction, per bin (Bech et al. 2003, eq. 2).

    y = terrain height − beam-center height at that range, same units as `a`.
    Returns 0 where y <= -a (beam fully clear), 1 where y >= a (beam fully
    blocked), and the circular-segment fraction in between. y == 0 gives
    exactly 0.5 (beam center grazing the terrain).
    """
    y = np.asarray(y, dtype=float)
    a = np.asarray(a, dtype=float)
    pbb = np.empty(np.broadcast(y, a).shape, dtype=float)

    clear = y <= -a
    block = y >= a
    partial = ~clear & ~block

    pbb[clear] = 0.0
    pbb[block] = 1.0

    yp, ap = np.broadcast_arrays(y, a)
    yp, ap = yp[partial], ap[partial]
    pbb[partial] = (
        yp * np.sqrt(ap**2 - yp**2) + ap**2 * np.arcsin(yp / ap) + np.pi * ap**2 / 2.0
    ) / (np.pi * ap**2)

    return pbb


def cum_beam_block_frac(pbb, axis=-1):
    """CBB: running-maximum PBB along range. See module docstring — never skip this."""
    return np.maximum.accumulate(pbb, axis=axis)
