from qgis.core import QgsProcessingProvider

from .alg_discover import DiscoverAlgorithm
from .alg_evaluate import EvaluateAlgorithm
from .alg_network import NetworkAlgorithm
from .alg_robustness import RobustnessAlgorithm


class RadarSitingToolkitProvider(QgsProcessingProvider):
    def id(self):
        return "radar_siting_toolkit"

    def name(self):
        return "LAMULA™ Radar Siting Toolkit"

    def loadAlgorithms(self):
        self.addAlgorithm(EvaluateAlgorithm())
        self.addAlgorithm(DiscoverAlgorithm())
        self.addAlgorithm(NetworkAlgorithm())
        self.addAlgorithm(RobustnessAlgorithm())
