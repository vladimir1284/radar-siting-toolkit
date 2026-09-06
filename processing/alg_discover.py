"""Discover mode (Plan_Tecnico_Radar_Siting_Toolkit.md, sections 2.2, 8, 9 -- Fase 2).

Takes an eligible-terrain mask and a DEM, runs a full 360-degree engine
sweep (core.engine.sweep_site) from every eligible mask cell with one
uniform assumed tower height, and writes a merit raster on the mask's own
grid: for each eligible cell, how much of the swept area is visible below
two height thresholds, plus its widest contiguous blocked sector.

This does NOT decide a site. Per section 2.2, it narrows a whole territory
down to a defensible short list for alg_evaluate.

Merit here is AREA-weighted (each polar sample weighted by its own range,
s*ds*dtheta, since equal azimuth/range steps cover more physical area at
longer range), not IMPORTANCE-weighted -- there is no population/importance
raster input here, deliberately: this mode's grid is a search space (every
eligible cell probed as a hypothetical site), not the area of interest that
an importance layer describes. The importance-weighted figure of merit
(section 4.1) is computed downstream, per short-listed candidate, in
alg_evaluate / alg_robustness.
"""

import os

import numpy as np
import rasterio
from qgis.core import (
    QgsProcessing,
    QgsProcessingAlgorithm,
    QgsProcessingException,
    QgsProcessingParameterNumber,
    QgsProcessingParameterRasterDestination,
    QgsProcessingParameterRasterLayer,
    QgsProcessingParameterString,
)
from qgis.PyQt.QtCore import QCoreApplication

from ..core.dem import Dem, DemNodataError, raster_cell_centers_lonlat
from ..core.engine import sweep_site
from ..core.manifest import build_manifest, write_manifest
from ..core.metrics import largest_contiguous_blocked_sector_deg, weighted_visible_fraction

NODATA_VALUE = -9999.0


class DiscoverAlgorithm(QgsProcessingAlgorithm):
    INPUT_DEM = "INPUT_DEM"
    INPUT_MASK = "INPUT_MASK"
    DEM_VERTICAL_DATUM = "DEM_VERTICAL_DATUM"
    TOWER_HEIGHT_M = "TOWER_HEIGHT_M"
    EARTH_RADIUS_FACTOR_K = "EARTH_RADIUS_FACTOR_K"
    ELEVATION_DEG = "ELEVATION_DEG"
    BEAMWIDTH_DEG = "BEAMWIDTH_DEG"
    MAX_RANGE_M = "MAX_RANGE_M"
    AZIMUTH_STEP_DEG = "AZIMUTH_STEP_DEG"
    VISIBLE_THRESHOLD_1_M = "VISIBLE_THRESHOLD_1_M"
    VISIBLE_THRESHOLD_2_M = "VISIBLE_THRESHOLD_2_M"
    BLOCKED_SECTOR_THRESHOLD = "BLOCKED_SECTOR_THRESHOLD"
    OUTPUT_MERIT = "OUTPUT_MERIT"

    def tr(self, string):
        return QCoreApplication.translate("DiscoverAlgorithm", string)

    def createInstance(self):
        return DiscoverAlgorithm()

    def name(self):
        return "discover"

    def displayName(self):
        return self.tr("Discover candidate sites")

    def group(self):
        return self.tr("Radar Siting Toolkit")

    def groupId(self):
        return "radarsitingtoolkit"

    def shortHelpString(self):
        return self.tr(
            "Runs the horizon-angle and Bech engine from every eligible cell of "
            "an eligible-terrain mask, with one assumed tower height for the "
            "whole search. Writes a 3-band merit raster on the mask's own grid: "
            "area-weighted visible fraction at two height thresholds, and widest "
            "contiguous blocked sector. Does not decide a site -- narrows the "
            "search space for alg_evaluate (section 2.2)."
        )

    def initAlgorithm(self, config=None):
        self.addParameter(QgsProcessingParameterRasterLayer(self.INPUT_DEM, self.tr("DEM")))
        self.addParameter(
            QgsProcessingParameterRasterLayer(
                self.INPUT_MASK, self.tr("Eligible-terrain mask (nonzero = eligible)")
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
            QgsProcessingParameterNumber(
                self.TOWER_HEIGHT_M,
                self.tr("Assumed tower height for the whole search [m]"),
                type=QgsProcessingParameterNumber.Double,
                defaultValue=20.0,
                minValue=0.0,
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
                self.tr("Azimuth step [deg] -- coarser is faster, less accurate"),
                type=QgsProcessingParameterNumber.Double,
                defaultValue=5.0,
                minValue=0.01,
                maxValue=45.0,
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
                self.tr("Height threshold 2 [m above local terrain]"),
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
            QgsProcessingParameterRasterDestination(self.OUTPUT_MERIT, self.tr("Merit raster"))
        )

    def processAlgorithm(self, parameters, context, feedback):
        dem_layer = self.parameterAsRasterLayer(parameters, self.INPUT_DEM, context)
        if dem_layer is None:
            raise QgsProcessingException(self.invalidRasterError(parameters, self.INPUT_DEM))
        mask_layer = self.parameterAsRasterLayer(parameters, self.INPUT_MASK, context)
        if mask_layer is None:
            raise QgsProcessingException(self.invalidRasterError(parameters, self.INPUT_MASK))

        vertical_datum = self.parameterAsString(parameters, self.DEM_VERTICAL_DATUM, context)
        tower_height_m = self.parameterAsDouble(parameters, self.TOWER_HEIGHT_M, context)
        k = self.parameterAsDouble(parameters, self.EARTH_RADIUS_FACTOR_K, context)
        elevation_deg = self.parameterAsDouble(parameters, self.ELEVATION_DEG, context)
        beamwidth_deg = self.parameterAsDouble(parameters, self.BEAMWIDTH_DEG, context)
        max_range_m = self.parameterAsDouble(parameters, self.MAX_RANGE_M, context)
        azimuth_step_deg = self.parameterAsDouble(parameters, self.AZIMUTH_STEP_DEG, context)
        threshold_1 = self.parameterAsDouble(parameters, self.VISIBLE_THRESHOLD_1_M, context)
        threshold_2 = self.parameterAsDouble(parameters, self.VISIBLE_THRESHOLD_2_M, context)
        blocked_sector_threshold = self.parameterAsDouble(
            parameters, self.BLOCKED_SECTOR_THRESHOLD, context
        )

        output_path = self.parameterAsOutputLayer(parameters, self.OUTPUT_MERIT, context)

        mask_path = mask_layer.source()
        with rasterio.open(mask_path) as mask_ds:
            mask_arr = mask_ds.read(1)
            mask_nodata = mask_ds.nodata
            transform = mask_ds.transform
            crs = mask_ds.crs
            height, width = mask_arr.shape

        valid = ~np.isnan(mask_arr.astype(float))
        if mask_nodata is not None:
            valid &= ~np.isclose(mask_arr.astype(float), mask_nodata)
        eligible = valid & (mask_arr != 0)

        lons, lats = raster_cell_centers_lonlat(transform, crs, height, width)

        merit_1 = np.full((height, width), NODATA_VALUE, dtype="float32")
        merit_2 = np.full((height, width), NODATA_VALUE, dtype="float32")
        blocked_sector = np.full((height, width), NODATA_VALUE, dtype="float32")

        dem_path = dem_layer.source()
        dem = Dem(dem_path, vertical_datum=vertical_datum or None)

        rows, cols = np.nonzero(eligible)
        total = len(rows) or 1
        try:
            for i, (row, col) in enumerate(zip(rows, cols)):
                if feedback.isCanceled():
                    break

                lon, lat = float(lons[row, col]), float(lats[row, col])
                try:
                    ground = dem.sample_lonlat([lon], [lat], nodata_policy="raise")[0]
                except DemNodataError:
                    feedback.pushInfo(f"cell ({row},{col}): no DEM data, left as nodata.")
                    continue

                h0 = ground + tower_height_m
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
                # Area weight, not importance weight -- see module docstring.
                area_weight = np.broadcast_to(out["s_m"], out["h_min_m"].shape)

                merit_1[row, col] = weighted_visible_fraction(out["h_min_m"], area_weight, threshold_1)
                merit_2[row, col] = weighted_visible_fraction(out["h_min_m"], area_weight, threshold_2)
                blocked_sector[row, col] = largest_contiguous_blocked_sector_deg(
                    out["cbb"][:, -1], azimuth_step_deg, threshold=blocked_sector_threshold
                )

                if i % 50 == 0:
                    feedback.setProgress(int(100 * i / total))
        finally:
            dem.close()

        profile = dict(
            driver="GTiff",
            height=height,
            width=width,
            count=3,
            dtype="float32",
            crs=crs,
            transform=transform,
            nodata=NODATA_VALUE,
        )
        with rasterio.open(output_path, "w", **profile) as dst:
            dst.write(merit_1, 1)
            dst.write(merit_2, 2)
            dst.write(blocked_sector, 3)
            dst.set_band_description(1, f"area_weighted_visible_fraction_at_{threshold_1:.0f}m")
            dst.set_band_description(2, f"area_weighted_visible_fraction_at_{threshold_2:.0f}m")
            dst.set_band_description(3, "largest_blocked_sector_deg")

        manifest = build_manifest(
            params={
                "tower_height_m": tower_height_m,
                "earth_radius_factor_k": k,
                "elevation_deg": elevation_deg,
                "beamwidth_deg": beamwidth_deg,
                "max_range_m": max_range_m,
                "azimuth_step_deg": azimuth_step_deg,
                "visible_threshold_1_m": threshold_1,
                "visible_threshold_2_m": threshold_2,
                "blocked_sector_threshold": blocked_sector_threshold,
                "dem_vertical_datum": vertical_datum or None,
            },
            input_layers={"dem": dem_path, "mask": mask_path},
        )
        manifest_path = os.path.splitext(output_path)[0] + "_manifest.json"
        write_manifest(manifest_path, manifest)

        feedback.pushInfo(f"Manifest written to {manifest_path}")

        return {self.OUTPUT_MERIT: output_path}
