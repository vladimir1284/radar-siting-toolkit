"""One-radar sweep: horizon angle, minimum visible height, and cumulative
beam blockage over a full 360-degree scan.

Plan_Tecnico_Radar_Siting_Toolkit.md, sections 3.1-3.3, 9 (Fase 1).

Factors out the per-azimuth loop that both the Fase 0 notebook and
processing/alg_evaluate.py need: sample the DEM along each ray at native
resolution, accumulate the horizon angle and cumulative beam-block
fraction, and stack the results into one (n_azimuths, n_ranges) grid per
quantity, for a single radar site.
"""

import numpy as np

from .blockage import beam_block_frac, beam_height, cum_beam_block_frac, half_power_radius
from .geometry import native_step, ray_points, transformer_from_aeqd, transformer_to_aeqd
from .horizon import horizon_angle, min_visible_height
from .importance import sample_importance_along_rays
from .metrics import largest_contiguous_blocked_sector_deg, weighted_visible_fraction


def sweep_site(
    dem,
    lon0,
    lat0,
    h0,
    k,
    elevation_deg,
    beamwidth_deg,
    max_range_m,
    azimuth_step_deg=1.0,
    range_step_m=None,
    nodata_policy="block",
):
    """Full 360-degree sweep from one radar site.

    dem : an open core.dem.Dem.
    lon0, lat0 : radar antenna position (WGS84).
    h0 : radar antenna height [m] = ground elevation at the radar + tower height.
    k : effective earth-radius factor.
    elevation_deg, beamwidth_deg : lowest scan elevation and beamwidth [deg].
    max_range_m : sweep range [m].
    azimuth_step_deg : azimuth sampling step [deg].
    range_step_m : DEM sampling step [m]; defaults to the DEM's own native
        resolution (section 3.5 -- never a coarser fixed step).
    nodata_policy : passed through to Dem.sample_lonlat ("block" by default,
        so a nodata cell blocks the beam instead of raising and aborting the
        whole sweep -- see core.dem).

    Returns a dict with azimuths_deg (n_azimuths,), s_m (n_ranges, shared by
    every ray), theta_rad / h_min_m / cbb / lons_deg / lats_deg, each
    (n_azimuths, n_ranges). lons_deg/lats_deg are each ray's own sample
    points in WGS84 -- e.g. for core.engine.site_merit to sample an
    importance layer at exactly the points this sweep already visited,
    without recomputing the ray geometry.
    """
    if range_step_m is None:
        range_step_m = native_step(dem.transform)

    to_lonlat = transformer_from_aeqd(lon0, lat0)
    azimuths_deg = np.arange(0.0, 360.0, azimuth_step_deg)

    elevation_rad = np.deg2rad(elevation_deg)
    beamwidth_rad = np.deg2rad(beamwidth_deg)

    s_ref = None
    theta_rows, hmin_rows, cbb_rows, lon_rows, lat_rows = [], [], [], [], []
    for az in azimuths_deg:
        s, x, y = ray_points(az, max_range_m, range_step_m)
        s_ref = s
        lons, lats = to_lonlat.transform(x, y)
        z = dem.sample_lonlat(lons, lats, nodata_policy=nodata_policy)

        theta = horizon_angle(z, s, h0, k)
        hmin = min_visible_height(z, s, h0, k, theta=theta)

        bh = beam_height(s, h0, elevation_rad, k)
        a = half_power_radius(s, beamwidth_rad)
        pbb = beam_block_frac(z - bh, a)
        cbb = cum_beam_block_frac(pbb)

        theta_rows.append(theta)
        hmin_rows.append(hmin)
        cbb_rows.append(cbb)
        lon_rows.append(lons)
        lat_rows.append(lats)

    return {
        "azimuths_deg": azimuths_deg,
        "s_m": s_ref,
        "theta_rad": np.array(theta_rows),
        "h_min_m": np.array(hmin_rows),
        "cbb": np.array(cbb_rows),
        "lons_deg": np.array(lon_rows),
        "lats_deg": np.array(lat_rows),
    }


def h_min_at_point(
    dem,
    radar_lon,
    radar_lat,
    h0,
    k,
    query_lon,
    query_lat,
    range_step_m=None,
    nodata_policy="block",
):
    """Exact minimum visible height [m] at one point, seen from one radar.

    Plan_Tecnico_Radar_Siting_Toolkit.md, section 2.3 (Fase 2, alg_network):
    joining several radars' coverage over an arbitrary area-of-interest grid
    needs h_min at points that rarely land on one of sweep_site's fixed
    360-azimuth rays. Rather than interpolate between two nearby rays, this
    samples the single ray from the radar to (query_lon, query_lat) at
    native DEM resolution, along the exact bearing to that point -- exact,
    at the cost of one full ray per (radar, query point) pair.

    Returns h0 itself if the query point is closer than one range step (no
    intervening terrain to accumulate a horizon over).
    """
    if range_step_m is None:
        range_step_m = native_step(dem.transform)

    to_xy = transformer_to_aeqd(radar_lon, radar_lat)
    x_q, y_q = to_xy.transform(query_lon, query_lat)
    s_total = float(np.hypot(x_q, y_q))
    if s_total < range_step_m:
        return h0

    az_deg = float(np.degrees(np.arctan2(x_q, y_q)) % 360.0)
    s, x, y = ray_points(az_deg, s_total, range_step_m)

    to_lonlat = transformer_from_aeqd(radar_lon, radar_lat)
    lons, lats = to_lonlat.transform(x, y)
    z = dem.sample_lonlat(lons, lats, nodata_policy=nodata_policy)

    theta = horizon_angle(z, s, h0, k)
    hmin = min_visible_height(z, s, h0, k, theta=theta)
    return float(hmin[-1])


def site_merit(
    dem,
    importance_band,
    importance_transform,
    importance_crs,
    importance_nodata,
    lon0,
    lat0,
    h0,
    k,
    elevation_deg,
    beamwidth_deg,
    max_range_m,
    azimuth_step_deg=1.0,
    thresholds=(1000.0, 3000.0),
    blocked_sector_threshold=0.5,
):
    """Importance-weighted figure of merit for one candidate site (section 4.1),
    plus its widest contiguous blocked sector (section 4.2 penalty) --
    the two headline numbers alg_robustness (Fase 3) sweeps over k, tower
    height and DEM source.

    Runs one full sweep_site scan and samples the importance layer at each
    ray's own points (core.importance.sample_importance_along_rays), then
    weights each sample by both the importance value there and the physical
    area it represents (proportional to range s, since equal azimuth/range
    steps cover more area at longer range -- same area correction
    alg_discover uses for its unweighted proxy, section 2.2). This needs no
    Cartesian resampling of the importance layer: the importance layer is
    read at the sweep's own points, not the other way around.

    Raises ValueError (core.metrics._weighted_fraction, propagated through
    weighted_visible_fraction) if this candidate's sweep does not overlap
    the importance layer at all -- e.g. a candidate sited far outside the
    area of interest. That is a real modeling error for this candidate/combo,
    not a bug to swallow; callers sweeping many candidates should catch it
    per candidate rather than aborting the whole matrix.
    """
    out = sweep_site(
        dem, lon0, lat0, h0, k, elevation_deg, beamwidth_deg, max_range_m, azimuth_step_deg
    )
    importance_at_rays = sample_importance_along_rays(
        importance_band,
        importance_transform,
        importance_crs,
        importance_nodata,
        out["lons_deg"],
        out["lats_deg"],
    )
    weight = importance_at_rays * out["s_m"][np.newaxis, :]

    merit = {
        threshold: weighted_visible_fraction(out["h_min_m"], weight, threshold)
        for threshold in thresholds
    }
    largest_blocked_sector_deg = largest_contiguous_blocked_sector_deg(
        out["cbb"][:, -1], azimuth_step_deg, threshold=blocked_sector_threshold
    )
    return {"merit": merit, "largest_blocked_sector_deg": largest_blocked_sector_deg}
