"""Evaluate mode (Plan_Tecnico_Radar_Siting_Toolkit.md, sections 2.1, 8, 9 -- Fase 1).

Takes a point layer of candidate sites (or existing radars) and a DEM, runs
the horizon + Bech engine (core.engine.sweep_site) over a full 360-degree
scan from each site, and writes:

- a copy of the input points with per-site summary fields attached;
- one per-azimuth CSV per site (core.metrics operates on the same
  per-azimuth summary, e.g. largest_contiguous_blocked_sector_deg);
- one execution manifest (core.manifest) for the whole run.

Importance-weighted merit scores (section 4.1, core.metrics.merit_scores)
need an area-of-interest raster grid, not just the per-site ray sweep, and
are out of scope here -- see alg_robustness (Fase 3 of the plan).
"""

import csv
import math
import os

import numpy as np
from qgis.core import (
    QgsCoordinateReferenceSystem,
    QgsCoordinateTransform,
    QgsFeature,
    QgsFeatureSink,
    QgsField,
    QgsFields,
    QgsGeometry,
    QgsPointXY,
    QgsProcessing,
    QgsProcessingAlgorithm,
    QgsProcessingException,
    QgsProcessingParameterFeatureSink,
    QgsProcessingParameterFeatureSource,
    QgsProcessingParameterFolderDestination,
    QgsProcessingParameterNumber,
    QgsProcessingParameterRasterLayer,
    QgsProcessingParameterString,
    QgsWkbTypes,
)
from qgis.PyQt.QtCore import QCoreApplication, QVariant

from ..core.dem import Dem, DemNodataError
from ..core.engine import sweep_site
from ..core.manifest import build_manifest, write_manifest
from ..core.metrics import largest_contiguous_blocked_sector_deg

WGS84 = QgsCoordinateReferenceSystem("EPSG:4326")


class EvaluateAlgorithm(QgsProcessingAlgorithm):
    INPUT_SITES = "INPUT_SITES"
    INPUT_DEM = "INPUT_DEM"
    DEM_VERTICAL_DATUM = "DEM_VERTICAL_DATUM"
    TOWER_HEIGHT_FIELD = "TOWER_HEIGHT_FIELD"
    GROUND_ELEVATION_FIELD = "GROUND_ELEVATION_FIELD"
    EARTH_RADIUS_FACTOR_K = "EARTH_RADIUS_FACTOR_K"
    ELEVATION_DEG = "ELEVATION_DEG"
    BEAMWIDTH_DEG = "BEAMWIDTH_DEG"
    MAX_RANGE_M = "MAX_RANGE_M"
    AZIMUTH_STEP_DEG = "AZIMUTH_STEP_DEG"
    BLOCKED_SECTOR_THRESHOLD = "BLOCKED_SECTOR_THRESHOLD"
    OUTPUT_SITES = "OUTPUT_SITES"
    OUTPUT_RAYS = "OUTPUT_RAYS"
    OUTPUT_FOLDER = "OUTPUT_FOLDER"

    def tr(self, string):
        return QCoreApplication.translate("EvaluateAlgorithm", string)

    def createInstance(self):
        return EvaluateAlgorithm()

    def name(self):
        return "evaluate"

    def displayName(self):
        return self.tr("Evaluate candidate sites")

    def group(self):
        return self.tr("Radar Siting Toolkit")

    def groupId(self):
        return "radarsitingtoolkit"

    def shortHelpString(self):
        return self.tr(
            "Runs the horizon-angle and Bech partial-beam-blockage engine over a "
            "full 360-degree scan from each candidate site (or existing radar), "
            "against a single DEM. Writes per-site summary fields, one "
            "per-azimuth CSV per site, and an execution manifest.\n\n"
            "Not included here: importance-weighted merit scores (section 4.1 of "
            "the plan) need an area-of-interest raster grid, not just the "
            "per-site ray sweep -- see the robustness-matrix algorithm."
        )

    def initAlgorithm(self, config=None):
        self.addParameter(
            QgsProcessingParameterFeatureSource(
                self.INPUT_SITES,
                self.tr("Candidate sites (or existing radars)"),
                [QgsProcessing.TypeVectorPoint],
            )
        )
        self.addParameter(
            QgsProcessingParameterRasterLayer(self.INPUT_DEM, self.tr("DEM"))
        )
        self.addParameter(
            QgsProcessingParameterString(
                self.DEM_VERTICAL_DATUM,
                self.tr("DEM vertical datum (recorded in the manifest only)"),
                defaultValue="",
                optional=True,
            )
        )
        self.addParameter(
            QgsProcessingParameterString(
                self.TOWER_HEIGHT_FIELD,
                self.tr("Field: tower height above ground [m]"),
                defaultValue="tower_height_m",
            )
        )
        self.addParameter(
            QgsProcessingParameterString(
                self.GROUND_ELEVATION_FIELD,
                self.tr(
                    "Field: ground elevation at the site [m] (optional -- "
                    "sampled from the DEM if absent or null)"
                ),
                defaultValue="ground_elevation_m",
                optional=True,
            )
        )
        self.addParameter(
            QgsProcessingParameterNumber(
                self.EARTH_RADIUS_FACTOR_K,
                self.tr("Effective earth-radius factor k"),
                type=QgsProcessingParameterNumber.Double,
                defaultValue=4.0 / 3.0,
                minValue=0.1,
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
                defaultValue=1.0,
                minValue=0.01,
                maxValue=45.0,
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
            QgsProcessingParameterFeatureSink(
                self.OUTPUT_SITES, self.tr("Site summary")
            )
        )
        self.addParameter(
            QgsProcessingParameterFeatureSink(
                self.OUTPUT_RAYS,
                self.tr(
                    "Per-azimuth blockage rays (site to max range, style by "
                    "cbb_at_max_range or min_clear_elev_deg)"
                ),
                type=QgsProcessing.TypeVectorLine,
                optional=True,
                createByDefault=True,
            )
        )
        self.addParameter(
            QgsProcessingParameterFolderDestination(
                self.OUTPUT_FOLDER, self.tr("Output folder (per-site CSV + manifest)")
            )
        )

    def processAlgorithm(self, parameters, context, feedback):
        source = self.parameterAsSource(parameters, self.INPUT_SITES, context)
        if source is None:
            raise QgsProcessingException(self.invalidSourceError(parameters, self.INPUT_SITES))

        dem_layer = self.parameterAsRasterLayer(parameters, self.INPUT_DEM, context)
        if dem_layer is None:
            raise QgsProcessingException(self.invalidRasterError(parameters, self.INPUT_DEM))
        vertical_datum = self.parameterAsString(parameters, self.DEM_VERTICAL_DATUM, context)

        tower_height_field = self.parameterAsString(parameters, self.TOWER_HEIGHT_FIELD, context)
        ground_elevation_field = self.parameterAsString(
            parameters, self.GROUND_ELEVATION_FIELD, context
        )
        k = self.parameterAsDouble(parameters, self.EARTH_RADIUS_FACTOR_K, context)
        elevation_deg = self.parameterAsDouble(parameters, self.ELEVATION_DEG, context)
        beamwidth_deg = self.parameterAsDouble(parameters, self.BEAMWIDTH_DEG, context)
        max_range_m = self.parameterAsDouble(parameters, self.MAX_RANGE_M, context)
        azimuth_step_deg = self.parameterAsDouble(parameters, self.AZIMUTH_STEP_DEG, context)
        blocked_sector_threshold = self.parameterAsDouble(
            parameters, self.BLOCKED_SECTOR_THRESHOLD, context
        )

        output_folder = self.parameterAsString(parameters, self.OUTPUT_FOLDER, context)
        os.makedirs(output_folder, exist_ok=True)

        fields = QgsFields(source.fields())
        for name, qvariant_type in (
            ("ground_elevation_m", QVariant.Double),
            ("h0_m", QVariant.Double),
            ("largest_blocked_sector_deg", QVariant.Double),
            ("worst_az_min_clear_elev_deg", QVariant.Double),
            ("worst_az_cbb_at_max_range", QVariant.Double),
            ("azimuth_csv", QVariant.String),
        ):
            fields.append(QgsField(name, qvariant_type))

        sink, dest_id = self.parameterAsSink(
            parameters,
            self.OUTPUT_SITES,
            context,
            fields,
            QgsWkbTypes.Point,
            source.sourceCrs(),
        )
        if sink is None:
            raise QgsProcessingException(self.invalidSinkError(parameters, self.OUTPUT_SITES))

        ray_fields = QgsFields()
        for name, qvariant_type in (
            ("site_name", QVariant.String),
            ("azimuth_deg", QVariant.Double),
            ("min_clear_elev_deg", QVariant.Double),
            ("h_min_m", QVariant.Double),
            ("cbb_at_max_range", QVariant.Double),
        ):
            ray_fields.append(QgsField(name, qvariant_type))
        ray_sink, ray_dest_id = self.parameterAsSink(
            parameters,
            self.OUTPUT_RAYS,
            context,
            ray_fields,
            QgsWkbTypes.LineString,
            source.sourceCrs(),
        )

        to_wgs84 = QgsCoordinateTransform(source.sourceCrs(), WGS84, context.transformContext())
        from_wgs84 = QgsCoordinateTransform(WGS84, source.sourceCrs(), context.transformContext())

        dem_path = dem_layer.source()
        dem = Dem(dem_path, vertical_datum=vertical_datum or None)

        site_params = {
            "earth_radius_factor_k": k,
            "elevation_deg": elevation_deg,
            "beamwidth_deg": beamwidth_deg,
            "max_range_m": max_range_m,
            "azimuth_step_deg": azimuth_step_deg,
            "blocked_sector_threshold": blocked_sector_threshold,
            "dem_vertical_datum": vertical_datum or None,
        }

        total = source.featureCount() or 1
        try:
            for current, feature in enumerate(source.getFeatures()):
                if feedback.isCanceled():
                    break

                geom = feature.geometry()
                if geom.isEmpty():
                    feedback.reportError(f"Feature {feature.id()}: empty geometry, skipped.")
                    continue
                point = geom.asPoint()
                point_wgs84 = to_wgs84.transform(point)
                lon, lat = point_wgs84.x(), point_wgs84.y()

                if tower_height_field not in feature.fields().names():
                    raise QgsProcessingException(
                        f"Field '{tower_height_field}' not found on the input layer."
                    )
                tower_height = feature[tower_height_field]
                if tower_height is None:
                    raise QgsProcessingException(
                        f"Feature {feature.id()}: '{tower_height_field}' is null."
                    )
                tower_height = float(tower_height)

                ground_elevation = None
                if ground_elevation_field and ground_elevation_field in feature.fields().names():
                    value = feature[ground_elevation_field]
                    if value is not None:
                        ground_elevation = float(value)
                if ground_elevation is None:
                    ground_elevation = float(
                        dem.sample_lonlat([lon], [lat], nodata_policy="raise")[0]
                    )

                h0 = ground_elevation + tower_height

                site_label = feature.attribute("name") if "name" in feature.fields().names() else None
                site_label = site_label or f"site_{feature.id()}"

                feedback.pushInfo(f"[{current + 1}/{total}] {site_label}: h0={h0:.1f} m")

                out = sweep_site(
                    dem,
                    lon0=lon,
                    lat0=lat,
                    h0=h0,
                    k=k,
                    elevation_deg=elevation_deg,
                    beamwidth_deg=beamwidth_deg,
                    max_range_m=max_range_m,
                    azimuth_step_deg=azimuth_step_deg,
                )

                cbb_edge = out["cbb"][:, -1]
                largest_sector = largest_contiguous_blocked_sector_deg(
                    cbb_edge, azimuth_step_deg, threshold=blocked_sector_threshold
                )

                if ray_sink is not None:
                    for i, az in enumerate(out["azimuths_deg"]):
                        end_lon = float(out["lons_deg"][i, -1])
                        end_lat = float(out["lats_deg"][i, -1])
                        ray_geom = QgsGeometry.fromPolylineXY(
                            [QgsPointXY(lon, lat), QgsPointXY(end_lon, end_lat)]
                        )
                        ray_geom.transform(from_wgs84)
                        ray_feature = QgsFeature(ray_fields)
                        ray_feature.setGeometry(ray_geom)
                        ray_feature.setAttributes(
                            [
                                site_label,
                                float(az),
                                math.degrees(float(out["theta_rad"][i, -1])),
                                float(out["h_min_m"][i, -1]),
                                float(out["cbb"][i, -1]),
                            ]
                        )
                        ray_sink.addFeature(ray_feature, QgsFeatureSink.FastInsert)

                azimuth_csv_path = os.path.join(output_folder, f"{_safe_filename(site_label)}.csv")
                _write_azimuth_summary_csv(
                    azimuth_csv_path,
                    azimuths_deg=out["azimuths_deg"],
                    theta_rad=out["theta_rad"],
                    h_min_m=out["h_min_m"],
                    cbb=out["cbb"],
                )

                out_feature = QgsFeature(fields)
                out_feature.setGeometry(geom)
                attrs = list(feature.attributes())
                # "Worst azimuth" = the direction with the highest minimum-clear-
                # elevation requirement (equivalently, the highest CBB) at max
                # range -- a one-number summary of the site's hardest-blocked
                # direction, not an average across azimuths.
                worst_idx = int(np.argmax(out["theta_rad"][:, -1]))
                attrs += [
                    ground_elevation,
                    h0,
                    largest_sector,
                    float(np.degrees(out["theta_rad"][worst_idx, -1])),
                    float(cbb_edge[worst_idx]),
                    azimuth_csv_path,
                ]
                out_feature.setAttributes(attrs)
                sink.addFeature(out_feature, QgsFeatureSink.FastInsert)

                feedback.setProgress(int(100 * (current + 1) / total))
        except DemNodataError as e:
            raise QgsProcessingException(str(e))
        finally:
            dem.close()

        manifest = build_manifest(
            params=site_params,
            input_layers={"dem": dem_path},
        )
        write_manifest(os.path.join(output_folder, "manifest.json"), manifest)

        result = {self.OUTPUT_SITES: dest_id, self.OUTPUT_FOLDER: output_folder}
        if ray_dest_id is not None:
            result[self.OUTPUT_RAYS] = ray_dest_id
        return result


def _safe_filename(label):
    return "".join(c if c.isalnum() or c in "-_" else "_" for c in str(label))


def _write_azimuth_summary_csv(path, azimuths_deg, theta_rad, h_min_m, cbb):
    """One row per azimuth: min clear elevation, h_min and CBB at max range.

    The full range-resolved grid is not written here (Fase 1 scope): see
    alg_robustness (Fase 3 of the plan) for CSV export driven by the
    robustness sweep.
    """
    with open(path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(
            [
                "azimuth_deg",
                "min_clear_elevation_deg_at_max_range",
                "h_min_m_at_max_range",
                "cbb_at_max_range",
            ]
        )
        for i, az in enumerate(azimuths_deg):
            writer.writerow(
                [
                    float(az),
                    math.degrees(float(theta_rad[i, -1])),
                    float(h_min_m[i, -1]),
                    float(cbb[i, -1]),
                ]
            )
