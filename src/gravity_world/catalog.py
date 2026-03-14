from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from .settings import Settings


@dataclass(frozen=True)
class IndicatorSpec:
    code: str
    slug: str
    date_from: int
    date_to: int


@dataclass(frozen=True)
class SourceSpec:
    id: str
    kind: str
    description: str
    output_subdir: str
    url: str | None = None
    filename: str | None = None
    landing_url: str | None = None
    base_url: str | None = None
    country_metadata_url: str | None = None
    format: str | None = None
    per_page: int | None = None
    expected_filenames: list[str] = field(default_factory=list)
    notes: str | None = None
    indicators: list[IndicatorSpec] = field(default_factory=list)


def _merge_dicts(base: dict[str, Any], override: dict[str, Any]) -> dict[str, Any]:
    merged = dict(base)
    for key, value in override.items():
        if isinstance(value, dict) and isinstance(merged.get(key), dict):
            merged[key] = _merge_dicts(merged[key], value)
        else:
            merged[key] = value
    return merged


def _index_by_id(items: list[dict[str, Any]]) -> dict[str, dict[str, Any]]:
    return {item["id"]: item for item in items}


def load_catalog(settings: Settings) -> list[SourceSpec]:
    catalog_path = settings.config_dir / "sources.json"
    overrides_path = settings.config_dir / "source_overrides.json"

    raw_catalog = json.loads(catalog_path.read_text(encoding="utf-8"))
    source_defs = raw_catalog["sources"]

    if overrides_path.exists():
        raw_overrides = json.loads(overrides_path.read_text(encoding="utf-8"))
        overrides_by_id = _index_by_id(raw_overrides.get("sources", []))
        merged_defs = []
        for source_def in source_defs:
            override = overrides_by_id.get(source_def["id"])
            merged_defs.append(_merge_dicts(source_def, override) if override else source_def)
        source_defs = merged_defs

    specs: list[SourceSpec] = []
    for source_def in source_defs:
        indicators = [
            IndicatorSpec(
                code=item["code"],
                slug=item["slug"],
                date_from=item["date_from"],
                date_to=item["date_to"],
            )
            for item in source_def.get("indicators", [])
        ]
        specs.append(
            SourceSpec(
                id=source_def["id"],
                kind=source_def["kind"],
                description=source_def["description"],
                output_subdir=source_def["output_subdir"],
                url=source_def.get("url"),
                filename=source_def.get("filename"),
                landing_url=source_def.get("landing_url"),
                base_url=source_def.get("base_url"),
                country_metadata_url=source_def.get("country_metadata_url"),
                format=source_def.get("format"),
                per_page=source_def.get("per_page"),
                expected_filenames=source_def.get("expected_filenames", []),
                notes=source_def.get("notes"),
                indicators=indicators,
            )
        )

    return specs


def get_source(settings: Settings, source_id: str) -> SourceSpec:
    for spec in load_catalog(settings):
        if spec.id == source_id:
            return spec
    raise KeyError(f"Unknown source id: {source_id}")
