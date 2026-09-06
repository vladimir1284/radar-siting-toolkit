"""Execution manifest: parameters, input checksums, versions, timestamp.

Plan_Tecnico_Radar_Siting_Toolkit.md, section 6.1.

Every run writes one of these next to its outputs. It exists to close a
single objection -- "how did you get this?" -- for a report that justifies
public investment: which inputs (by content, not just by path, since a path
can be silently overwritten), which parameters, which library versions, and
when.
"""

import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

_PACKAGE_ROOT = Path(__file__).resolve().parent.parent
_METADATA_PATH = _PACKAGE_ROOT / "metadata.txt"


def sha256_file(path, chunk_size=1 << 20):
    """Streamed SHA-256 hex digest of a file, safe for multi-GB DEMs."""
    import hashlib

    digest = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def plugin_version(metadata_path=_METADATA_PATH):
    """Plugin version string, read from metadata.txt's `version=` line."""
    with open(metadata_path, encoding="utf-8") as f:
        for line in f:
            if line.startswith("version="):
                return line.split("=", 1)[1].strip()
    raise ValueError(f"no version= line found in {metadata_path}")


def gdal_version():
    """GDAL version as seen through rasterio, or 'unknown' if rasterio is absent.

    Inside QGIS, GDAL is the application's own bundled GDAL; rasterio is a
    thin wrapper over the same library, so this reports the version actually
    used to read the DEM (core.dem), not some unrelated system GDAL.
    """
    try:
        import rasterio

        return rasterio.__gdal_version__
    except ImportError:
        return "unknown"


def build_manifest(params, input_layers, plugin_version_override=None):
    """Assemble the manifest dict for one run.

    params : JSON-serializable dict of every parameter that affects the
        result (k, tower height, DEM source choice, thresholds, ...).
        Serialization is deliberately left to fail loudly on a
        non-serializable value rather than silently stringifying it.
    input_layers : {name: path} for every input raster/vector layer
        (DEM(s), importance layer, candidate sites). Each is hashed, so a
        manifest proves *which* file version was used, not just which path.
    """
    layers = {}
    for name, path in input_layers.items():
        path = Path(path)
        layers[name] = {"path": str(path), "sha256": sha256_file(path)}

    return {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "plugin_version": plugin_version_override or plugin_version(),
        "numpy_version": np.__version__,
        "gdal_version": gdal_version(),
        "parameters": params,
        "input_layers": layers,
    }


def write_manifest(path, manifest):
    """Write a manifest dict as pretty-printed JSON to `path`."""
    Path(path).write_text(json.dumps(manifest, indent=2))
