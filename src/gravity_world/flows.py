from __future__ import annotations

from pathlib import Path

from .harmonize import harmonize_source_frame
from .settings import Settings


FLOW_METHOD_COLUMNS = [
    "da_min_open",
    "da_min_closed",
    "da_pb_closed",
    "da_pb_open",
    "ffr_closed",
    "ffr_open",
]

REQUIRED_ID_COLUMNS = ["year0", "orig", "dest"]


def _require_pandas():
    try:
        import pandas as pd  # type: ignore
    except ImportError as exc:  # pragma: no cover
        raise RuntimeError("pandas is required for flow normalization commands.") from exc
    return pd


def _build_code_map(frame, settings: Settings, raw_column: str, prefix: str):
    unique_codes = frame[[raw_column]].drop_duplicates().copy()
    unique_codes[raw_column] = unique_codes[raw_column].astype("string")

    harmonized, unmatched = harmonize_source_frame(
        unique_codes,
        settings,
        source_id="abel_cohen_flows",
        code_column=raw_column,
    )

    mapping = harmonized[
        [
            raw_column,
            "canonical_id",
            "iso3",
            "canonical_name",
            "entity_type",
            "is_current",
            "match_method",
        ]
    ].copy()
    mapping = mapping.rename(
        columns={
            "canonical_id": f"{prefix}_canonical_id",
            "iso3": f"{prefix}_iso3",
            "canonical_name": f"{prefix}_name",
            "entity_type": f"{prefix}_entity_type",
            "is_current": f"{prefix}_is_current",
            "match_method": f"{prefix}_match_method",
        }
    )

    counts = frame.groupby(raw_column, dropna=False).size().reset_index(name="row_count")
    unmatched = unmatched[[raw_column]].drop_duplicates().merge(counts, on=raw_column, how="left")
    unmatched = unmatched.rename(columns={raw_column: f"{prefix}_code_raw"})
    return mapping, unmatched


def normalize_abel_cohen_flows(
    settings: Settings,
    preferred_flow: str = "da_pb_closed",
    raw_path: Path | None = None,
    output_dir: Path | None = None,
) -> list[Path]:
    pd = _require_pandas()

    if preferred_flow not in FLOW_METHOD_COLUMNS:
        allowed = ", ".join(FLOW_METHOD_COLUMNS)
        raise ValueError(f"Unknown preferred flow column: {preferred_flow}. Allowed values: {allowed}")

    raw_path = raw_path or settings.raw_dir / "abel_cohen" / "abel_cohen_flows_1990_2020.csv"
    if not raw_path.exists():
        raise FileNotFoundError(
            "Abel-Cohen raw flow file is missing. Run `gravity-data download --source abel_cohen_flows` first."
        )
    if raw_path.stat().st_size == 0:
        raise ValueError(
            "The Abel-Cohen raw flow file is empty. Re-download it outside the sandbox or replace it manually before normalization."
        )

    frame = pd.read_csv(raw_path)
    missing_columns = [column for column in REQUIRED_ID_COLUMNS if column not in frame.columns]
    if missing_columns:
        actual = ", ".join(frame.columns)
        missing = ", ".join(missing_columns)
        raise ValueError(f"Abel-Cohen file is missing required columns: {missing}. Actual columns: {actual}")

    available_flow_columns = [column for column in FLOW_METHOD_COLUMNS if column in frame.columns]
    if not available_flow_columns:
        actual = ", ".join(frame.columns)
        expected = ", ".join(FLOW_METHOD_COLUMNS)
        raise ValueError(
            f"Abel-Cohen file does not contain any supported flow columns. Expected one of: {expected}. Actual columns: {actual}"
        )
    if preferred_flow not in available_flow_columns:
        available = ", ".join(available_flow_columns)
        raise ValueError(
            f"Preferred flow column `{preferred_flow}` is not present in the Abel-Cohen file. Available flow columns: {available}"
        )

    working = frame[REQUIRED_ID_COLUMNS + available_flow_columns].copy()
    working["source_dataset"] = "abel_cohen_flows"
    working["period_start_year"] = pd.to_numeric(working["year0"], errors="coerce").astype("Int64")
    working["period_end_year"] = working["period_start_year"] + 5
    working["period_label"] = working["period_start_year"].astype("string") + "-" + working["period_end_year"].astype("string")
    working = working.rename(columns={"orig": "origin_code_raw", "dest": "destination_code_raw"})

    for column in available_flow_columns:
        working[column] = pd.to_numeric(working[column], errors="coerce")
    for column in FLOW_METHOD_COLUMNS:
        if column not in working.columns:
            working[column] = pd.NA

    origin_map, origin_unmatched = _build_code_map(working, settings, "origin_code_raw", "origin")
    destination_map, destination_unmatched = _build_code_map(working, settings, "destination_code_raw", "destination")

    normalized = working.merge(origin_map, on="origin_code_raw", how="left")
    normalized = normalized.merge(destination_map, on="destination_code_raw", how="left")
    normalized["is_self_flow"] = normalized["origin_canonical_id"].eq(normalized["destination_canonical_id"])
    normalized["flow"] = normalized[preferred_flow]
    normalized["flow_measure"] = preferred_flow

    output_dir = output_dir or settings.processed_dir / "flows"
    output_dir.mkdir(parents=True, exist_ok=True)

    normalized_path = output_dir / "bilateral_flows_5y_abel_cohen.csv"
    origin_unmatched_path = output_dir / "abel_cohen_unmatched_origins.csv"
    destination_unmatched_path = output_dir / "abel_cohen_unmatched_destinations.csv"

    output_columns = [
        "source_dataset",
        "period_start_year",
        "period_end_year",
        "period_label",
        "origin_code_raw",
        "origin_canonical_id",
        "origin_iso3",
        "origin_name",
        "origin_entity_type",
        "origin_is_current",
        "origin_match_method",
        "destination_code_raw",
        "destination_canonical_id",
        "destination_iso3",
        "destination_name",
        "destination_entity_type",
        "destination_is_current",
        "destination_match_method",
        *FLOW_METHOD_COLUMNS,
        "flow",
        "flow_measure",
        "is_self_flow",
    ]
    normalized = normalized[output_columns].sort_values(
        ["period_start_year", "origin_code_raw", "destination_code_raw"]
    ).reset_index(drop=True)

    normalized.to_csv(normalized_path, index=False)
    origin_unmatched.to_csv(origin_unmatched_path, index=False)
    destination_unmatched.to_csv(destination_unmatched_path, index=False)
    return [normalized_path, origin_unmatched_path, destination_unmatched_path]