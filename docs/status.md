# Scope & Status

## Out of scope

Explicitly excluded, so expectations stay defensible in front of a
government or donor reviewer:

- Rain attenuation and frequency-band selection — a first-order decision in
  the tropics, and not resolvable with geometry alone.
- Real anomalous propagation (ducting, measured refractivity gradients) —
  only the effective earth-radius factor `k` is parameterized.
- Sea clutter, ground clutter, RF interference, spectrum compatibility.
- Quantitative precipitation estimation, nowcasting, radar data processing.
- Tower structure, radome, wind load, electrical resilience.
- Cost, lifecycle, staffing.

**This toolkit covers roughly 15% of a full siting study.**

## Status

Early scaffolding. Per the project plan, the analysis engine (`core/`) was
developed and validated in a notebook first — the critical path for the
originating Jamaica field report — with `core/` and `processing/` (the QGIS
plugin wrapper: Evaluate, Discover, Network, Robustness matrix) following
afterward. All four Processing algorithms are implemented; QGIS Plugin
Manager packaging/publication has not happened yet (`metadata.txt`:
`experimental=True`).

See `Plan_Tecnico_Radar_Siting_Toolkit.md` (repository root) for the full
technical plan, phase breakdown, and the reasoning behind each design
decision referenced throughout this documentation.

## Testing

`tests/` covers every `core/` module (`pytest`, configured via `pytest.ini`
with `pythonpath = .`). `notebooks/fase0_motor_horizonte_bech.ipynb`
validates the engine against `wradlib` on the real Jamaica radar site and
DEM data (Nivel 2 validation, section 7 of the plan) — `wradlib` is a
dev/notebook-only dependency (`requirements-dev.txt`), never a plugin
runtime dependency.

```sh
pip install -r requirements.txt -r requirements-dev.txt
pytest
```

## License

GPL v2 or later (`LICENSE`) — the de facto requirement for QGIS plugins.
