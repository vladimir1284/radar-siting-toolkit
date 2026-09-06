"""Network mode (Plan_Tecnico_Radar_Siting_Toolkit.md, sections 2.3, 4.2, 8, 9 -- Fase 2).

Takes a set of radars (own candidates plus existing radars, domestic or
foreign) and an importance layer, and computes joint coverage and residual
gaps over the importance layer's own grid: for every radar, the exact
minimum visible height at every importance-layer cell within range
(core.engine.h_min_at_point -- one exact ray per radar/cell pair, not an
interpolation between two of sweep_site's fixed azimuths, since cells
essentially never land on one), then combines them.

Unlike alg_discover, this uses the importance layer's real weights (section
4.3): the output is the importance-weighted figure of merit for a *given*
network configuration, not a search-space proxy.
"""

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
    QgsProcessingParameterFeatureSource,
    QgsProcessingParameterNumber,
    QgsProcessingParameterRasterDestination,
    QgsProcessingParameterRasterLayer,
    QgsProcessingParameterString,
)
from qgis.PyQt.QtCore import QCoreApplication

from ..core.dem import Dem, DemNodataError, raster_cell_centers_lonlat
from ..core.engine import h_min_at_point
from ..core.geometry import transformer_to_aeqd
from ..core.manifest import build_manifest, write_manifest
from ..core.metrics import network_residual_gap_fraction, weighted_visible_fraction

NODATA_VALUE = -9999.0
WGS84 = QgsCoordinateReferenceSystem("EPSG:4326")


class NetworkAlgorithm(QgsProcessingAlgorithm):
    INPUT_RADARS = "INPUT_RADARS"
    INPUT_DEM = "INPUT_DEM"
    INPUT_IMPORTANCE = "INPUT_IMPORTANCE"
    DEM_VERTICAL_DATUM = "DEM_VERTICAL_DATUM"
    TOWER_HEIGHT_FIELD = "TOWER_HEIGHT_FIELD"
    GROUND_ELEVATION_FIELD = "GROUND_ELEVATION_FIELD"
    EARTH_RADIUS_FACTOR_K = "EARTH_RADIUS_FACTOR_K"
    MAX_RANGE_M = "MAX_RANGE_M"
    VISIBLE_THRESHOLD_1_M = "VISIBLE_THRESHOLD_1_M"
    VISIBLE_THRESHOLD_2_M = "VISIBLE_THRESHOLD_2_M"
    OUTPUT_BEST_HMIN = "OUTPUT_BEST_HMIN"

    def tr(self, string):
        return QCoreApplication.translate("NetworkAlgorithm", string)

    def createInstance(self):
        return NetworkAlgorithm()

    def name(self):
        return "network"

    def displayName(self):
        return self.tr("Network coverage")

    def group(self):
        return self.tr("Radar Siting Toolkit")

    def groupId(self):
        return "radarsitingtoolkit"

    def shortHelpString(self):
        return self.tr(
            "Joint coverage across a set of radars (candidates plus existing "
            "radars), over the importance layer's own grid. For every radar and "
            "every importance-layer cell within range, computes the exact "
            "minimum visible height (core.engine.h_min_at_point), takes the "
            "per-cell minimum across radars, and reports importance-weighted "
            "joint visible fraction at two height thresholds and the residual "
            "gap fraction (area seen by no radar) at the second threshold "
            "(section 4.2). Writes the per-cell best (minimum) h_min as a "
            "raster, and a JSON summary alongside the manifest."
        )

    def initAlgorithm(self, config=None):
        self.addParameter(
            QgsProcessingParameterFeatureSource(
                self.INPUT_RADARS,
                self.tr("Radars (candidates and/or existing)"),
                [QgsProcessing.TypeVectorPoint],
            )
        )
        self.addParameter(QgsProcessingParameterRasterLayer(self.INPUT_DEM, self.tr("DEM")))
        self.addParameter(
            QgsProcessingParameterRasterLayer(
                self.INPUT_IMPORTANCE,
                self.tr("Importance layer (area of interest + weights, section 4.3)"),
            )
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
                self.MAX_RANGE_M,
                self.tr("Max range per radar [m]"),
                type=QgsProcessingParameterNumber.Double,
                defaultValue=250_000.0,
                minValue=1.0,
            )
        )
        self.addParameter(
            QgsProcessingParameterNumber(
                self.VISIBLE_THRESHOLD_1_M,
                self.tr("Height threshold 1 [m above local terrain]"),
                type=QgsProcessingParameterNumber.Double,
                defaultValue=1000.0,
                minValue=0.0,
            )
        )
        self.addParameter(
            QgsProcessingParameterNumber(
                self.VISIBLE_THRESHOLD_2_M,
                self.tr(
                    "Height threshold 2 [m above local terrain] -- also used for "
                    "the residual-gap fraction (section 4.2)"
                ),
                type=QgsProcessingParameterNumber.Double,
                defaultValue=3000.0,
                minValue=0.0,
            )
        )
        self.addParameter(
            QgsProcessingParameterRasterDestination(
                self.OUTPUT_BEST_HMIN, self.tr("Best (joint) minimum visible height")
            )
        )

    def processAlgorithm(self, parameters, context, feedback):
        source = self.parameterAsSource(parameters, self.INPUT_RADARS, context)
        if source is None:
            raise QgsProcessingException(self.invalidSourceError(parameters, self.INPUT_RADARS))
        dem_layer = self.parameterAsRasterLayer(parameters, self.INPUT_DEM, context)
        if dem_layer is None:
            raise QgsProcessingException(self.invalidRasterError(parameters, self.INPUT_DEM))
        importance_layer = self.parameterAsRasterLayer(parameters, self.INPUT_IMPORTANCE, context)
        if importance_layer is None:
            raise QgsProcessingException(self.invalidRasterError(parameters, self.INPUT_IMPORTANCE))

        vertical_datum = self.parameterAsString(parameters, self.DEM_VERTICAL_DATUM, context)
        tower_height_field = self.parameterAsString(parameters, self.TOWER_HEIGHT_FIELD, context)
        ground_elevation_field = self.parameterAsString(
            parameters, self.GROUND_ELEVATION_FIELD, context
        )
        k = self.parameterAsDouble(parameters, self.EARTH_RADIUS_FACTOR_K, context)
        max_range_m = self.parameterAsDouble(parameters, self.MAX_RANGE_M, context)
        threshold_1 = self.parameterAsDouble(parameters, self.VISIBLE_THRESHOLD_1_M, context)
        threshold_2 = self.parameterAsDouble(parameters, self.VISIBLE_THRESHOLD_2_M, context)

        output_path = self.parameterAsOutputLayer(parameters, self.OUTPUT_BEST_HMIN, context)

        importance_path = importance_layer.source()
        with rasterio.open(importance_path) as imp_ds:
            importance_arr = imp_ds.read(1).astype(float)
            imp_nodata = imp_ds.nodata
            transform = imp_ds.transform
            crs = imp_ds.crs
            height, width = importance_arr.shape

        valid = ~np.isnan(importance_arr)
        if imp_nodata is not None:
            valid &= ~np.isclose(importance_arr, imp_nodata)
        importance_arr = np.where(valid, importance_arr, 0.0)

        lons, lats = raster_cell_centers_lonlat(transform, crs, height, width)

        dem_path = dem_layer.source()
        dem = Dem(dem_path, vertical_datum=vertical_datum or None)

        to_wgs84 = QgsCoordinateTransform(source.sourceCrs(), WGS84, context.transformContext())

        radars = []
        for feature in source.getFeatures():
            geom = feature.geometry()
            if geom.isEmpty():
                feedback.reportError(f"Feature {feature.id()}: empty geometry, skipped.")
                continue
            point_wgs84 = to_wgs84.transform(geom.asPoint())
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
                ground_elevation = float(dem.sample_lonlat([lon], [lat], nodata_policy="raise")[0])

            label = feature.attribute("name") if "name" in feature.fields().names() else None
            label = label or f"radar_{feature.id()}"
            radars.append(
                {"label": label, "lon": lon, "lat": lat, "h0": ground_elevation + tower_height}
            )

        if not radars:
            raise QgsProcessingException("No radars found in the input layer.")

        h_min_stack = np.full((len(radars), height, width), np.inf, dtype=float)
        rows, cols = np.nonzero(valid)
        total = len(radars) * (len(rows) or 1)
        done = 0

        try:
            for r_idx, radar in enumerate(radars):
                feedback.pushInfo(f"Radar {r_idx + 1}/{len(radars)}: {radar['label']}")
                to_xy = transformer_to_aeqd(radar["lon"], radar["lat"])
                cell_lons = lons[rows, cols]
                cell_lats = lats[rows, cols]
                xs, ys = to_xy.transform(cell_lons, cell_lats)
                distances = np.hypot(xs, ys)

                for row, col, lon, lat, s_total in zip(rows, cols, cell_lons, cell_lats, distances):
                    if feedback.isCanceled():
                        break
                    done += 1
                    if done % 200 == 0:
                        feedback.setProgress(int(100 * done / total))
                    if s_total > max_range_m:
                        continue
                    try:
                        h_min_stack[r_idx, row, col] = h_min_at_point(
                            dem,
                            radar["lon"],
                            radar["lat"],
                            radar["h0"],
                            k,
                            float(lon),
                            float(lat),
                        )
                    except DemNodataError:
                        continue
                if feedback.isCanceled():
                    break
        finally:
            dem.close()

        best_h_min = np.min(h_min_stack, axis=0)

        joint_visible_1 = weighted_visible_fraction(best_h_min, importance_arr, threshold_1, mask=valid)
        joint_visible_2 = weighted_visible_fraction(best_h_min, importance_arr, threshold_2, mask=valid)
        residual_gap = network_residual_gap_fraction(
            h_min_stack, importance_arr, threshold_m=threshold_2, mask=valid, axis=0
        )

        feedback.pushInfo(
            f"Joint visible fraction: {threshold_1:.0f} m = {joint_visible_1:.3f}, "
            f"{threshold_2:.0f} m = {joint_visible_2:.3f}. "
            f"Residual gap fraction (> {threshold_2:.0f} m, no radar): {residual_gap:.3f}"
        )

        out_arr = np.where(valid & np.isfinite(best_h_min), best_h_min, NODATA_VALUE).astype("float32")
        profile = dict(
            driver="GTiff",
            height=height,
            width=width,
            count=1,
            dtype="float32",
            crs=crs,
            transform=transform,
            nodata=NODATA_VALUE,
        )
        with rasterio.open(output_path, "w", **profile) as dst:
            dst.write(out_arr, 1)
            dst.set_band_description(1, "best_h_min_m_across_radars")

        summary = {
            "radars": [{"label": r["label"], "lon": r["lon"], "lat": r["lat"], "h0_m": r["h0"]} for r in radars],
            "visible_threshold_1_m": threshold_1,
            "joint_visible_fraction_1": joint_visible_1,
            "visible_threshold_2_m": threshold_2,
            "joint_visible_fraction_2": joint_visible_2,
            "residual_gap_fraction_at_threshold_2": residual_gap,
        }
        summary_path = os.path.splitext(output_path)[0] + "_summary.json"
        with open(summary_path, "w") as f:
            json.dump(summary, f, indent=2)

        manifest = build_manifest(
            params={
                "earth_radius_factor_k": k,
                "max_range_m": max_range_m,
                "visible_threshold_1_m": threshold_1,
                "visible_threshold_2_m": threshold_2,
                "dem_vertical_datum": vertical_datum or None,
                "radars": [r["label"] for r in radars],
            },
            input_layers={"dem": dem_path, "importance": importance_path},
        )
        manifest_path = os.path.splitext(output_path)[0] + "_manifest.json"
        write_manifest(manifest_path, manifest)

        feedback.pushInfo(f"Summary written to {summary_path}, manifest to {manifest_path}")

        return {self.OUTPUT_BEST_HMIN: output_path}
