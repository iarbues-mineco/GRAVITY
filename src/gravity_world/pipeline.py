from __future__ import annotations

from pathlib import Path

from .catalog import SourceSpec, load_catalog
from .settings import Settings
from .sources.manual_file import ManualFileSource
from .sources.static_file import StaticFileSource
from .sources.world_bank import WorldBankSource


def inventory(settings: Settings) -> list[dict[str, str]]:
    records = []
    for spec in load_catalog(settings):
        automation = "automatic"
        if spec.kind == "manual_file" and not spec.url:
            automation = "manual_or_override"
        records.append(
            {
                "id": spec.id,
                "kind": spec.kind,
                "automation": automation,
                "description": spec.description,
            }
        )
    return records


def run_downloads(settings: Settings, source_ids: list[str] | None = None) -> dict[str, list[Path]]:
    settings.ensure_directories()
    selected = []
    for spec in load_catalog(settings):
        if source_ids and spec.id not in source_ids:
            continue
        selected.append(spec)

    outputs: dict[str, list[Path]] = {}
    for spec in selected:
        outputs[spec.id] = _run_source(settings, spec)
    return outputs


def _run_source(settings: Settings, spec: SourceSpec) -> list[Path]:
    if spec.kind == "world_bank":
        return WorldBankSource(spec, settings.raw_dir).download()
    if spec.kind == "static_file":
        return [StaticFileSource(spec, settings.raw_dir).download()]
    if spec.kind == "manual_file":
        return [ManualFileSource(spec, settings.raw_dir).prepare()]
    raise ValueError(f"Unsupported source kind: {spec.kind}")
