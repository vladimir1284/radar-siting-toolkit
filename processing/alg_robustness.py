"""Robustness matrix (Plan_Tecnico_Radar_Siting_Toolkit.md, sections 4, 6.2,
6.3, 9 -- Fase 3).

For every candidate site, sweeps the importance-weighted figure of merit
(core.engine.site_merit, section 4.1) across every combination of DEM
source, effective earth-radius factor k, and tower height, then ranks
candidates within each combination. If the ranking is identical across
every combination, the recommendation is shielded (section 6.2); if it
flips, the team finds out before a reviewer does.

Ranking is driven by the higher height threshold (threshold 2, e.g. 3 km,
typically the more demanding "can it see cyclone structure" number) -- both
thresholds are still reported per row so a reviewer can re-rank by the
other one by hand.

Also writes the marginal merit gained per additional meter of tower height
(section 6.2's budget argument) and, if matplotlib is importable, publication
figures (section 6.3). matplotlib is optional and lazily imported: this
algorithm's hard runtime dependencies stay numpy/rasterio/pyproj, same as
the rest of core/ (section 3.3) -- CSV/JSON output never depends on it.
"""

import csv
import json
import os

import numpy as np
import rasterio
from qgis.core import (
    QgsCoordinateReferenceSystem,
    QgsCoordinateTransform,
    QgsProcessing,
    QgsProcessingAlgorithm,
    QgsProcessingException,
    QgsProcessingParameterBoolean,
    QgsProcessingParameterFeatureSource,
    QgsProcessingParameterFolderDestination,
    QgsProcessingParameterMultipleLayers,
    QgsProcessingParameterNumber,
    QgsProcessingParameterRasterLayer,
    QgsProcessingParameterString,
)
from qgis.PyQt.QtCore import QCoreApplication

from ..core.dem import Dem, DemNodataError
from ..core.engine import site_merit
from ..core.manifest import build_manifest, write_manifest

WGS84 = QgsCoordinateReferenceSystem("EPSG:4326")


def _parse_float_list(text, param_name):
    try:
        values = [float(v.strip()) for v in text.split(",") if v.strip()]
    except ValueError:
        raise QgsProcessingException(f"{param_name}: could not parse '{text}' as a comma-separated list of numbers.")
    if not values:
        raise QgsProcessingException(f"{param_name}: at least one value is required.")
    return values


class RobustnessAlgorithm(QgsProcessingAlgorithm):
    INPUT_SITES = "INPUT_SITES"
    GROUND_ELEVATION_FIELD = "GROUND_ELEVATION_FIELD"
    INPUT_DEM_SOURCES = "INPUT_DEM_SOURCES"
    INPUT_IMPORTANCE = "INPUT_IMPORTANCE"
    K_VALUES = "K_VALUES"
    TOWER_HEIGHTS_M = "TOWER_HEIGHTS_M"
    ELEVATION_DEG = "ELEVATION_DEG"
    BEAMWIDTH_DEG = "BEAMWIDTH_DEG"
    MAX_RANGE_M = "MAX_RANGE_M"
    AZIMUTH_STEP_DEG = "AZIMUTH_STEP_DEG"
    VISIBLE_THRESHOLD_1_M = "VISIBLE_THRESHOLD_1_M"
    VISIBLE_THRESHOLD_2_M = "VISIBLE_THRESHOLD_2_M"
    BLOCKED_SECTOR_THRESHOLD = "BLOCKED_SECTOR_THRESHOLD"
    MAKE_FIGURES = "MAKE_FIGURES"
    OUTPUT_FOLDER = "OUTPUT_FOLDER"

    def tr(self, string):
        return QCoreApplication.translate("RobustnessAlgorithm", string)

    def createInstance(self):
        return RobustnessAlgorithm()

    def name(self):
        return "robustness"

    def displayName(self):
        return self.tr("Robustness matrix")

    def group(self):
        return self.tr("Radar Siting Toolkit")

    def groupId(self):
        return "radarsitingtoolkit"

    def shortHelpString(self):
        return self.tr(
            "Sweeps the importance-weighted figure of merit (section 4.1) for "
            "every candidate site across every combination of DEM source, "
            "effective earth-radius factor k, and tower height, and ranks "
            "candidates within each combination by the higher height threshold. "
            "Writes the full comparison table, the marginal merit gained per "
            "additional meter of tower height, a stability verdict, a manifest, "
            "and (if matplotlib is available) publication figures."
        )

    def initAlgorithm(self, config=None):
        self.addParameter(
            QgsProcessingParameterFeatureSource(
                self.INPUT_SITES, self.tr("Candidate sites"), [QgsProcessing.TypeVectorPoint]
            )
        )
        self.addParameter(
            QgsProcessingParameterString(
                self.GROUND_ELEVATION_FIELD,
                self.tr(
                    "Field: ground elevation at the site [m] (optional -- "
                    "sampled per DEM source if absent or null)"
                ),
                defaultValue="ground_elevation_m",
                optional=True,
            )
        )
        self.addParameter(
            QgsProcessingParameterMultipleLayers(
                self.INPUT_DEM_SOURCES,
                self.tr("DEM sources to sweep (e.g. GLO-30, FABDEM)"),
                layerType=QgsProcessing.TypeRaster,
            )
        )
        self.addParameter(
            QgsProcessingParameterRasterLayer(
                self.INPUT_IMPORTANCE,
                self.tr("Importance layer (area of interest + weights, section 4.3)"),
            )
        )
        self.addParameter(
            QgsProcessingParameterString(
                self.K_VALUES,
                self.tr("Effective earth-radius factors k to sweep (comma-separated)"),
                defaultValue="1.0,1.3333,2.0",
            )
        )
        self.addParameter(
            QgsProcessingParameterString(
                self.TOWER_HEIGHTS_M,
                self.tr("Tower heights to sweep [m] (comma-separated)"),
                defaultValue="10,15,20,25,30",
            )
        )
        self.addParameter(
            QgsProcessingParameterNumber(
                self.ELEVATION_DEG,
                self.tr("Lowest scan elevation [deg]"),
                type=QgsProcessingParameterNumber.Double,
                defaultValue=0.5,
            )
        )
        self.addParameter(
            QgsProcessingParameterNumber(
                self.BEAMWIDTH_DEG,
                self.tr("Beamwidth [deg]"),
                type=QgsProcessingParameterNumber.Double,
                defaultValue=1.0,
                minValue=0.01,
            )
        )
        self.addParameter(
            QgsProcessingParameterNumber(
                self.MAX_RANGE_M,
                self.tr("Max range [m]"),
                type=QgsProcessingParameterNumber.Double,
                defaultValue=250_000.0,
                minValue=1.0,
            )
        )
        self.addParameter(
            QgsProcessingParameterNumber(
                self.AZIMUTH_STEP_DEG,
                self.tr("Azimuth step [deg]"),
                type=QgsProcessingParameterNumber.Double,
                defaultValue=2.0,
                minValue=0.01,
                maxValue=45.0,
            )
        )
        self.addParameter(
            QgsProcessingParameterNumber(
                self.VISIBLE_THRESHOLD_1_M,
                self.tr("Height threshold 1 [m above local terrain] -- e.g. flash-flood QPE"),
                type=QgsProcessingParameterNumber.Double,
                defaultValue=1000.0,
                minValue=0.0,
            )
        )
        self.addParameter(
            QgsProcessingParameterNumber(
                self.VISIBLE_THRESHOLD_2_M,
                self.tr(
                    "Height threshold 2 [m above local terrain] -- e.g. cyclone "
                    "structure; drives the ranking"
                ),
                type=QgsProcessingParameterNumber.Double,
                defaultValue=3000.0,
                minValue=0.0,
            )
        )
        self.addParameter(
            QgsProcessingParameterNumber(
                self.BLOCKED_SECTOR_THRESHOLD,
                self.tr("CBB threshold counted as 'blocked' for the widest-sector metric"),
                type=QgsProcessingParameterNumber.Double,
                defaultValue=0.5,
                minValue=0.0,
                maxValue=1.0,
            )
        )
        self.addParameter(
            QgsProcessingParameterBoolean(
                self.MAKE_FIGURES,
                self.tr("Write publication figures (requires matplotlib)"),
                defaultValue=True,
            )
        )
        self.addParameter(
            QgsProcessingParameterFolderDestination(
                self.OUTPUT_FOLDER, self.tr("Output folder (CSV + JSON + manifest [+ figures])")
            )
        )

    def processAlgorithm(self, parameters, context, feedback):
        source = self.parameterAsSource(parameters, self.INPUT_SITES, context)
        if source is None:
            raise QgsProcessingException(self.invalidSourceError(parameters, self.INPUT_SITES))
        ground_elevation_field = self.parameterAsString(
            parameters, self.GROUND_ELEVATION_FIELD, context
        )
        dem_layers = self.parameterAsLayerList(parameters, self.INPUT_DEM_SOURCES, context)
        if not dem_layers:
            raise QgsProcessingException("At least one DEM source is required.")
        importance_layer = self.parameterAsRasterLayer(parameters, self.INPUT_IMPORTANCE, context)
        if importance_layer is None:
            raise QgsProcessingException(self.invalidRasterError(parameters, self.INPUT_IMPORTANCE))

        k_values = _parse_float_list(
            self.parameterAsString(parameters, self.K_VALUES, context), "k values"
        )
        tower_heights = _parse_float_list(
            self.parameterAsString(parameters, self.TOWER_HEIGHTS_M, context), "tower heights"
        )
        elevation_deg = self.parameterAsDouble(parameters, self.ELEVATION_DEG, context)
        beamwidth_deg = self.parameterAsDouble(parameters, self.BEAMWIDTH_DEG, context)
        max_range_m = self.parameterAsDouble(parameters, self.MAX_RANGE_M, context)
        azimuth_step_deg = self.parameterAsDouble(parameters, self.AZIMUTH_STEP_DEG, context)
        threshold_1 = self.parameterAsDouble(parameters, self.VISIBLE_THRESHOLD_1_M, context)
        threshold_2 = self.parameterAsDouble(parameters, self.VISIBLE_THRESHOLD_2_M, context)
        blocked_sector_threshold = self.parameterAsDouble(
            parameters, self.BLOCKED_SECTOR_THRESHOLD, context
        )
        make_figures = self.parameterAsBoolean(parameters, self.MAKE_FIGURES, context)

        output_folder = self.parameterAsString(parameters, self.OUTPUT_FOLDER, context)
        os.makedirs(output_folder, exist_ok=True)

        importance_path = importance_layer.source()
        with rasterio.open(importance_path) as imp_ds:
            importance_band = imp_ds.read(1).astype(float)
            importance_transform = imp_ds.transform
            importance_crs = imp_ds.crs
            importance_nodata = imp_ds.nodata

        to_wgs84 = QgsCoordinateTransform(source.sourceCrs(), WGS84, context.transformContext())
        sites = []
        for feature in source.getFeatures():
            geom = feature.geometry()
            if geom.isEmpty():
                feedback.reportError(f"Feature {feature.id()}: empty geometry, skipped.")
                continue
            point_wgs84 = to_wgs84.transform(geom.asPoint())
            label = feature.attribute("name") if "name" in feature.fields().names() else None
            label = label or f"site_{feature.id()}"

            fixed_ground_elevation = None
            if ground_elevation_field and ground_elevation_field in feature.fields().names():
                value = feature[ground_elevation_field]
                if value is not None:
                    fixed_ground_elevation = float(value)

            sites.append(
                {
                    "label": label,
                    "lon": point_wgs84.x(),
                    "lat": point_wgs84.y(),
                    "fixed_ground_elevation_m": fixed_ground_elevation,
                }
            )
        if not sites:
            raise QgsProcessingException("No candidate sites found in the input layer.")

        rows = []
        total = len(dem_layers) * len(sites) * len(k_values) * len(tower_heights) or 1
        done = 0

        for dem_layer in dem_layers:
            dem_label = dem_layer.name()
            dem = Dem(dem_layer.source())
            try:
                for site in sites:
                    if feedback.isCanceled():
                        break
                    ground_elevation = site["fixed_ground_elevation_m"]
                    if ground_elevation is None:
                        try:
                            ground_elevation = float(
                                dem.sample_lonlat([site["lon"]], [site["lat"]], nodata_policy="raise")[0]
                            )
                        except DemNodataError:
                            feedback.reportError(
                                f"{site['label']}/{dem_label}: no DEM data at the site, skipped."
                            )
                            done += len(k_values) * len(tower_heights)
                            continue

                    for k in k_values:
                        for tower_height in tower_heights:
                            if feedback.isCanceled():
                                break
                            done += 1
                            if done % 10 == 0:
                                feedback.setProgress(int(100 * done / total))

                            h0 = ground_elevation + tower_height
                            try:
                                result = site_merit(
                                    dem,
                                    importance_band,
                                    importance_transform,
                                    importance_crs,
                                    importance_nodata,
                                    site["lon"],
                                    site["lat"],
                                    h0,
                                    k,
                                    elevation_deg,
                                    beamwidth_deg,
                                    max_range_m,
                                    azimuth_step_deg,
                                    thresholds=(threshold_1, threshold_2),
                                    blocked_sector_threshold=blocked_sector_threshold,
                                )
                            except ValueError as e:
                                feedback.reportError(
                                    f"{site['label']}/{dem_label}/k={k}/tower={tower_height}: {e}"
                                )
                                merit_1 = merit_2 = blocked_sector = float("nan")
                            else:
                                merit_1 = result["merit"][threshold_1]
                                merit_2 = result["merit"][threshold_2]
                                blocked_sector = result["largest_blocked_sector_deg"]

                            rows.append(
                                {
                                    "site": site["label"],
                                    "dem_source": dem_label,
                                    "k": k,
                                    "tower_height_m": tower_height,
                                    "h0_m": h0,
                                    "merit_threshold_1": merit_1,
                                    "merit_threshold_2": merit_2,
                                    "largest_blocked_sector_deg": blocked_sector,
                                }
                            )
            finally:
                dem.close()

        _rank_within_combinations(rows)
        stability = _rank_stability(rows)
        marginal = _marginal_performance(rows)

        matrix_path = os.path.join(output_folder, "robustness_matrix.csv")
        _write_csv(
            matrix_path,
            rows,
            [
                "site", "dem_source", "k", "tower_height_m", "h0_m",
                "merit_threshold_1", "merit_threshold_2",
                "largest_blocked_sector_deg", "rank",
            ],
        )

        marginal_path = os.path.join(output_folder, "marginal_performance.csv")
        _write_csv(
            marginal_path,
            marginal,
            [
                "site", "dem_source", "k", "tower_height_from_m", "tower_height_to_m",
                "delta_merit_threshold_2_per_m",
            ],
        )

        summary = {
            "n_sites": len(sites),
            "n_combinations": len(dem_layers) * len(k_values) * len(tower_heights),
            "ranking_threshold_m": threshold_2,
            "ranking_stable_across_all_combinations": stability["order_stable"],
            "single_top_site_across_all_combinations": stability["top_site"],
            "rank_range_by_site": stability["rank_range_by_site"],
        }
        summary_path = os.path.join(output_folder, "summary.json")
        with open(summary_path, "w") as f:
            json.dump(summary, f, indent=2)

        feedback.pushInfo(
            "Ranking stable across all combinations: "
            f"{stability['order_stable']}. Top site: {stability['top_site'] or 'varies'}."
        )

        if make_figures:
            try:
                _write_figures(output_folder, rows, stability, threshold_2)
            except ImportError:
                feedback.pushInfo(
                    "matplotlib not available -- skipped figures. "
                    "CSV/JSON outputs are unaffected."
                )

        manifest = build_manifest(
            params={
                "k_values": k_values,
                "tower_heights_m": tower_heights,
                "elevation_deg": elevation_deg,
                "beamwidth_deg": beamwidth_deg,
                "max_range_m": max_range_m,
                "azimuth_step_deg": azimuth_step_deg,
                "visible_threshold_1_m": threshold_1,
                "visible_threshold_2_m": threshold_2,
                "blocked_sector_threshold": blocked_sector_threshold,
                "dem_sources": [layer.name() for layer in dem_layers],
                "sites": [s["label"] for s in sites],
            },
            input_layers={
                "importance": importance_path,
                **{f"dem_{layer.name()}": layer.source() for layer in dem_layers},
            },
        )
        write_manifest(os.path.join(output_folder, "manifest.json"), manifest)

        return {self.OUTPUT_FOLDER: output_folder}


def _rank_within_combinations(rows):
    """Assign `rank` in place: 1 = best merit_threshold_2 within the same
    (dem_source, k, tower_height_m) combination. NaN merit (a candidate whose
    sweep failed for that combination) always ranks last."""
    by_combo = {}
    for row in rows:
        key = (row["dem_source"], row["k"], row["tower_height_m"])
        by_combo.setdefault(key, []).append(row)

    for combo_rows in by_combo.values():
        ordered = sorted(
            combo_rows,
            key=lambda r: (np.isnan(r["merit_threshold_2"]), -_nan_to_neg_inf(r["merit_threshold_2"])),
        )
        for rank, row in enumerate(ordered, start=1):
            row["rank"] = rank


def _nan_to_neg_inf(value):
    return float("-inf") if np.isnan(value) else value


def _rank_stability(rows):
    orders_by_combo = {}
    top_by_combo = {}
    rank_range_by_site = {}

    for row in rows:
        combo_key = (row["dem_source"], row["k"], row["tower_height_m"])
        orders_by_combo.setdefault(combo_key, []).append((row["rank"], row["site"]))
        if row["rank"] == 1:
            top_by_combo[combo_key] = row["site"]
        prev = rank_range_by_site.get(row["site"], (row["rank"], row["rank"]))
        rank_range_by_site[row["site"]] = (min(prev[0], row["rank"]), max(prev[1], row["rank"]))

    orders = set()
    for combo_key, ranked in orders_by_combo.items():
        order = tuple(site for _, site in sorted(ranked))
        orders.add(order)

    top_sites = set(top_by_combo.values())

    return {
        "order_stable": len(orders) <= 1,
        "top_site": next(iter(top_sites)) if len(top_sites) == 1 else None,
        "rank_range_by_site": {site: list(r) for site, r in rank_range_by_site.items()},
    }


def _marginal_performance(rows):
    by_line = {}
    for row in rows:
        key = (row["site"], row["dem_source"], row["k"])
        by_line.setdefault(key, []).append(row)

    marginal = []
    for (site, dem_source, k), line_rows in by_line.items():
        line_rows = sorted(line_rows, key=lambda r: r["tower_height_m"])
        for prev_row, next_row in zip(line_rows, line_rows[1:]):
            delta_tower = next_row["tower_height_m"] - prev_row["tower_height_m"]
            if delta_tower == 0:
                continue
            delta_merit = next_row["merit_threshold_2"] - prev_row["merit_threshold_2"]
            marginal.append(
                {
                    "site": site,
                    "dem_source": dem_source,
                    "k": k,
                    "tower_height_from_m": prev_row["tower_height_m"],
                    "tower_height_to_m": next_row["tower_height_m"],
                    "delta_merit_threshold_2_per_m": delta_merit / delta_tower,
                }
            )
    return marginal


def _write_csv(path, rows, fieldnames):
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def _write_figures(output_folder, rows, stability, threshold_2):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figures_dir = os.path.join(output_folder, "figures")
    os.makedirs(figures_dir, exist_ok=True)

    sites = sorted({row["site"] for row in rows})
    for site in sites:
        fig, ax = plt.subplots(figsize=(7, 4), dpi=200)
        site_rows = [r for r in rows if r["site"] == site]
        by_line = {}
        for row in site_rows:
            key = (row["dem_source"], row["k"])
            by_line.setdefault(key, []).append(row)
        for (dem_source, k), line_rows in sorted(by_line.items()):
            line_rows = sorted(line_rows, key=lambda r: r["tower_height_m"])
            ax.plot(
                [r["tower_height_m"] for r in line_rows],
                [r["merit_threshold_2"] for r in line_rows],
                marker="o",
                label=f"{dem_source}, k={k:.3g}",
            )
        ax.set_xlabel("Tower height [m]")
        ax.set_ylabel(f"Importance-weighted visible fraction (H <= {threshold_2:.0f} m)")
        ax.set_title(f"Merit vs. tower height -- {site}")
        ax.legend(fontsize="small")
        ax.set_ylim(0, 1.05)
        fig.tight_layout()
        fig.savefig(os.path.join(figures_dir, f"merit_vs_tower_{_safe_filename(site)}.png"))
        plt.close(fig)

    fig, ax = plt.subplots(figsize=(7, 4), dpi=200)
    rank_ranges = stability["rank_range_by_site"]
    labels = sorted(rank_ranges.keys())
    lows = [rank_ranges[s][0] for s in labels]
    highs = [rank_ranges[s][1] for s in labels]
    ax.bar(labels, [h - l + 1 for l, h in zip(lows, highs)], bottom=lows)
    ax.set_ylabel("Rank across all combinations (1 = best)")
    ax.set_title("Rank stability by site")
    ax.invert_yaxis()
    fig.tight_layout()
    fig.savefig(os.path.join(figures_dir, "rank_stability.png"))
    plt.close(fig)


def _safe_filename(label):
    return "".join(c if c.isalnum() or c in "-_" else "_" for c in str(label))
