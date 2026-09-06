"""Merit figure and penalties: importance-weighted coverage, blockage geometry.

Plan_Tecnico_Radar_Siting_Toolkit.md, section 4.

Every function here is a plain array reduction with no I/O: callers are
responsible for producing `h_min` (core.horizon.min_visible_height) and
`importance` on a shared grid over the area of interest, and for supplying
an explicit `mask` when some cells are not part of that area (or are
nodata-blocked per core.dem). There is no default that silently treats
missing weight as zero importance for the whole grid -- an all-zero or
fully masked-out weight raster raises, per the same fail-loud policy as
core.dem, rather than returning a meaningless 0/0 fraction.
"""

from itertools import groupby

import numpy as np


def _weighted_fraction(indicator, importance, mask=None):
    """Importance-weighted mean of a boolean/float indicator over valid cells."""
    indicator = np.asarray(indicator, dtype=float)
    importance = np.asarray(importance, dtype=float)
    if indicator.shape != importance.shape:
        raise ValueError(
            f"indicator shape {indicator.shape} != importance shape {importance.shape}"
        )
    if mask is None:
        mask = np.ones(indicator.shape, dtype=bool)
    else:
        mask = np.asarray(mask, dtype=bool)
        if mask.shape != indicator.shape:
            raise ValueError(f"mask shape {mask.shape} != indicator shape {indicator.shape}")

    weight = importance * mask
    total_weight = weight.sum()
    if not np.isfinite(total_weight) or total_weight <= 0:
        raise ValueError(
            "importance layer has zero (or non-finite) total weight over the "
            "valid cells -- nothing to weight the fraction by. Section 4.3: "
            "the importance layer is a mandatory input, not a placeholder."
        )
    return float((indicator * weight).sum() / total_weight)


def weighted_visible_fraction(h_min, importance, threshold_m, mask=None):
    """Fraction of importance-weighted area where h_min <= threshold_m (section 4.1).

    h_min, importance : arrays of identical shape over the area of interest.
    threshold_m : height threshold [m above local terrain], e.g. 1000 or 3000.
    mask : optional boolean array, True for cells that count (inside the area
        of interest and not nodata-blocked). Defaults to all cells.
    """
    h_min = np.asarray(h_min, dtype=float)
    visible = h_min <= threshold_m
    return _weighted_fraction(visible, importance, mask)


def merit_scores(h_min, importance, thresholds=(1000.0, 3000.0), mask=None):
    """The two (or more) headline numbers per candidate (section 4.1).

    Returns {threshold_m: weighted_visible_fraction} -- deliberately never
    collapsed into one score: a site can be excellent for cyclone structure
    (3 km) and poor for flash-flood QPE (1 km), or vice versa.
    """
    return {
        threshold: weighted_visible_fraction(h_min, importance, threshold, mask)
        for threshold in thresholds
    }


def _max_circular_run(blocked):
    """Length of the longest run of True in a circular boolean sequence."""
    n = len(blocked)
    if n == 0 or not np.any(blocked):
        return 0
    if np.all(blocked):
        return n
    doubled = np.concatenate([blocked, blocked])
    longest = max(len(list(g)) for k, g in groupby(doubled) if k)
    return min(longest, n)


def largest_contiguous_blocked_sector_deg(cbb_by_azimuth, azimuth_step_deg, threshold=0.5):
    """Widest contiguous azimuthal sector with CBB > threshold (section 4.2).

    cbb_by_azimuth : (n_azimuths,) array, one cumulative-blockage value per
        azimuth (e.g. CBB at the outer edge of the area of interest -- CBB is
        a running max along range, so its value at the last range in scope
        already summarizes "blocked by then" for that azimuth). Azimuths are
        treated as wrapping circularly (0 deg adjacent to 360 deg): a 20%
        blockage scattered in pinches is benign, the same 20% concentrated
        in one 70 deg sector is a permanent hole, and a sector straddling
        due north must not be split into two by an arbitrary array boundary.
    """
    cbb_by_azimuth = np.asarray(cbb_by_azimuth, dtype=float)
    blocked = cbb_by_azimuth > threshold
    return _max_circular_run(blocked) * azimuth_step_deg


def network_residual_gap_fraction(h_min_stack, importance, threshold_m=3000.0, mask=None, axis=0):
    """Importance-weighted area seen by NO radar below threshold_m (section 4.2).

    h_min_stack : (n_radars, ...) array of per-radar h_min stacked along
        `axis`, all sharing the same grid as `importance`. A cell is a
        residual gap when every radar's h_min there exceeds threshold_m --
        the investment does not resolve it regardless of how many radars
        are in the configuration.
    """
    h_min_stack = np.asarray(h_min_stack, dtype=float)
    best_h_min = np.min(h_min_stack, axis=axis)
    gap = best_h_min > threshold_m
    return _weighted_fraction(gap, importance, mask)
