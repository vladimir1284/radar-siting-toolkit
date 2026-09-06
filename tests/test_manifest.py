import json

import numpy as np
import pytest

from core.manifest import build_manifest, plugin_version, sha256_file, write_manifest


def test_sha256_file_matches_hashlib(tmp_path):
    import hashlib

    p = tmp_path / "dem.tif"
    p.write_bytes(b"not really a geotiff, just some bytes" * 1000)
    expected = hashlib.sha256(p.read_bytes()).hexdigest()
    assert sha256_file(p) == expected


def test_plugin_version_reads_metadata_txt():
    # metadata.txt version=0.0.1 at time of writing; assert it's non-empty
    # and matches the file directly, not a hardcoded string that would go
    # stale the next time the plugin version bumps.
    from pathlib import Path

    metadata_path = Path(__file__).resolve().parent.parent / "metadata.txt"
    expected = None
    for line in metadata_path.read_text().splitlines():
        if line.startswith("version="):
            expected = line.split("=", 1)[1].strip()
    assert expected is not None
    assert plugin_version() == expected


def test_plugin_version_raises_without_version_line(tmp_path):
    bad_metadata = tmp_path / "metadata.txt"
    bad_metadata.write_text("[general]\nname=x\n")
    with pytest.raises(ValueError):
        plugin_version(metadata_path=bad_metadata)


def test_build_manifest_hashes_inputs_and_records_versions(tmp_path):
    dem_path = tmp_path / "dem.tif"
    dem_path.write_bytes(b"fake dem bytes")
    importance_path = tmp_path / "importance.tif"
    importance_path.write_bytes(b"fake importance bytes")

    manifest = build_manifest(
        params={"k": 4 / 3, "tower_height_m": 20},
        input_layers={"dem": dem_path, "importance": importance_path},
        plugin_version_override="9.9.9",
    )

    assert manifest["plugin_version"] == "9.9.9"
    assert manifest["numpy_version"] == np.__version__
    assert manifest["parameters"] == {"k": 4 / 3, "tower_height_m": 20}
    assert manifest["input_layers"]["dem"]["sha256"] == sha256_file(dem_path)
    assert manifest["input_layers"]["importance"]["sha256"] == sha256_file(importance_path)
    assert "timestamp_utc" in manifest


def test_write_manifest_round_trips_through_json(tmp_path):
    manifest = {"a": 1, "b": [1, 2, 3]}
    out = tmp_path / "manifest.json"
    write_manifest(out, manifest)
    assert json.loads(out.read_text()) == manifest


def test_build_manifest_rejects_non_serializable_params(tmp_path):
    dem_path = tmp_path / "dem.tif"
    dem_path.write_bytes(b"x")
    manifest = build_manifest(
        params={"bad": object()},
        input_layers={"dem": dem_path},
        plugin_version_override="9.9.9",
    )
    out = tmp_path / "manifest.json"
    with pytest.raises(TypeError):
        write_manifest(out, manifest)
