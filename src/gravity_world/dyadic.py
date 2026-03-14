from __future__ import annotations

from pathlib import Path
from zipfile import ZipFile

from .harmonize import harmonize_source_frame
from .settings import Settings


REQUIRED_CEPII_COLUMNS = [
    "iso_o",
    "iso_d",
    "contig",
    "comlang_off",
    "comlang_ethno",
    "colony",
    "comcol",
    "curcol",
    "col45",
    "smctry",
    "dist",
    "distcap",
    "distw",
    "distwces",
]


def _require_pandas():
    try:
        import pandas as pd  # type: ignore
    except ImportError as exc:  # pragma: no cover
        raise RuntimeError("pandas is required for CEPII normalization commands.") from exc
    return pd


def _ensure_cepii_xls(settings: Settings) -> Path:
    zip_path = settings.raw_dir / "cepii" / "dist_cepii.zip"
    if not zip_path.exists():
        raise FileNotFoundError(
            "CEPII GeoDist archive is missing. Place `dist_cepii.zip` in `data/raw/cepii/` first."
        )

    extract_dir = settings.interim_dir / "cepii" / "unzipped"
    extract_dir.mkdir(parents=True, exist_ok=True)
    xls_path = extract_dir / "dist_cepii.xls"
    if not xls_path.exists():
        with ZipFile(zip_path) as archive:
            archive.extract("dist_cepii.xls", path=extract_dir)
    return xls_path


def _build_code_map(frame, settings: Settings, raw_column: str, prefix: str):
    unique_codes = frame[[raw_column]].drop_duplicates().copy()
    unique_codes[raw_column] = unique_codes[raw_column].astype("string")

    harmonized, unmatched = harmonize_source_frame(
        unique_codes,
        settings,
        source_id="cepii_gravity",
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
    return mapping, unmatched.rename(columns={raw_column: f"{prefix}_code_raw"})


def normalize_cepii_controls(settings: Settings, output_dir: Path | None = None) -> list[Path]:
    pd = _require_pandas()
    import numpy as np

    xls_path = _ensure_cepii_xls(settings)
    frame = pd.read_excel(xls_path)

    missing = [column for column in REQUIRED_CEPII_COLUMNS if column not in frame.columns]
    if missing:
        actual = ", ".join(frame.columns)
        needed = ", ".join(missing)
        raise ValueError(f"CEPII file is missing required columns: {needed}. Actual columns: {actual}")

    working = frame[REQUIRED_CEPII_COLUMNS].copy()
    working["iso_o"] = working["iso_o"].astype("string").str.strip().str.upper()
    working["iso_d"] = working["iso_d"].astype("string").str.strip().str.upper()

    origin_map, origin_unmatched = _build_code_map(working, settings, "iso_o", "origin")
    destination_map, destination_unmatched = _build_code_map(working, settings, "iso_d", "destination")

    normalized = working.merge(origin_map, on="iso_o", how="left")
    normalized = normalized.merge(destination_map, on="iso_d", how="left")

    for column in ["dist", "distcap", "distw", "distwces"]:
        normalized[column] = pd.to_numeric(normalized[column], errors="coerce")
    for column in ["contig", "comlang_off", "comlang_ethno", "colony", "comcol", "curcol", "col45", "smctry"]:
        normalized[column] = pd.to_numeric(normalized[column], errors="coerce").fillna(0).astype(int)

    normalized["distance_km"] = normalized["distw"]
    normalized["distance_measure"] = "distw"
    normalized["ln_distance_km"] = np.where(normalized["distance_km"] > 0, np.log(normalized["distance_km"]), np.nan)

    output_dir = output_dir or settings.processed_dir / "dyadic"
    output_dir.mkdir(parents=True, exist_ok=True)

    controls_path = output_dir / "cepii_geodist_controls.csv"
    origin_unmatched_path = output_dir / "cepii_unmatched_origins.csv"
    destination_unmatched_path = output_dir / "cepii_unmatched_destinations.csv"

    output_columns = [
        "origin_canonical_id",
        "origin_iso3",
        "origin_name",
        "origin_entity_type",
        "origin_is_current",
        "origin_match_method",
        "destination_canonical_id",
        "destination_iso3",
        "destination_name",
        "destination_entity_type",
        "destination_is_current",
        "destination_match_method",
        "distance_measure",
        "distance_km",
        "ln_distance_km",
        "contig",
        "comlang_off",
        "comlang_ethno",
        "colony",
        "comcol",
        "curcol",
        "col45",
        "smctry",
        "dist",
        "distcap",
        "distw",
        "distwces",
    ]
    normalized = normalized[output_columns].sort_values(["origin_canonical_id", "destination_canonical_id"]).reset_index(drop=True)
    normalized.to_csv(controls_path, index=False)
    origin_unmatched.drop_duplicates().to_csv(origin_unmatched_path, index=False)
    destination_unmatched.drop_duplicates().to_csv(destination_unmatched_path, index=False)
    return [controls_path, origin_unmatched_path, destination_unmatched_path]