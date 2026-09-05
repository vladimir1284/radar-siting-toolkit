from qgis.core import QgsApplication

from .processing.provider import RadarSitingToolkitProvider


class RadarSitingToolkitPlugin:
    def __init__(self):
        self.provider = None

    def initGui(self):
        self.provider = RadarSitingToolkitProvider()
        QgsApplication.processingRegistry().addProvider(self.provider)

    def unload(self):
        QgsApplication.processingRegistry().removeProvider(self.provider)
