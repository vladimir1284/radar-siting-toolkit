from qgis.core import QgsProcessingProvider


class RadarSitingToolkitProvider(QgsProcessingProvider):
    def id(self):
        return "radar_siting_toolkit"

    def name(self):
        return "Radar Siting Toolkit"

    def loadAlgorithms(self):
        # Phase 1-3: register alg_evaluate, alg_discover, alg_network, alg_robustness here.
        pass
