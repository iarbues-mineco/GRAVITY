from __future__ import annotations

from pathlib import Path

from .settings import Settings


REFERENCE_COLUMNS = [
    "canonical_id",
    "canonical_name",
    "iso3",
    "iso2",
    "entity_type",
    "is_current",
    "valid_from",
    "valid_to",
    "world_bank_name",
    "world_bank_region",
    "world_bank_income_level",
    "notes",
]

OVERRIDE_COLUMNS = [
    "source_id",
    "source_code",
    "source_name",
    "canonical_id",
    "canonical_name",
    "iso3",
    "iso2",
    "entity_type",
    "is_current",
    "match_priority",
    "notes",
]


def _require_pandas():
    try:
        import pandas as pd  # type: ignore
    except ImportError as exc:  # pragma: no cover
        raise RuntimeError("pandas is required for harmonization commands.") from exc
    return pd


def _normalize_series(series):
    return (
        series.astype("string")
        .str.strip()
        .str.replace(r"\s+", " ", regex=True)
        .str.upper()
    )


def _load_world_bank_metadata(settings: Settings):
    pd = _require_pandas()
    path = settings.raw_dir / "world_bank" / "country_metadata.csv"
    if not path.exists():
        raise FileNotFoundError(
            "World Bank country metadata is missing. Run `gravity-data download --source world_bank_wdi` first."
        )
    return pd.read_csv(path, dtype="string")


def _load_manual_countries(settings: Settings):
    pd = _require_pandas()
    return pd.read_csv(settings.config_dir / "harmonization" / "manual_countries.csv", dtype="string")


def _load_source_overrides(settings: Settings):
    pd = _require_pandas()
    return pd.read_csv(settings.config_dir / "harmonization" / "source_overrides.csv", dtype="string")


def build_country_reference(settings: Settings) -> list[Path]:
    pd = _require_pandas()
    settings.ensure_directories()

    metadata = _load_world_bank_metadata(settings)
    manual = _load_manual_countries(settings)
    overrides = _load_source_overrides(settings)

    base = metadata.loc[metadata["region_name"].fillna("").ne("Aggregates")].copy()
    base["canonical_id"] = _normalize_series(base["id"])
    base["canonical_name"] = base["name"].astype("string").str.strip()
    base["iso3"] = base["canonical_id"]
    base["iso2"] = base["iso2c"].astype("string").str.strip().str.upper()
    base["entity_type"] = "country"
    base["is_current"] = "1"
    base["valid_from"] = pd.NA
    base["valid_to"] = pd.NA
    base["world_bank_name"] = base["name"].astype("string").str.strip()
    base["world_bank_region"] = base["region_name"].astype("string").str.strip()
    base["world_bank_income_level"] = base["income_level_name"].astype("string").str.strip()
    base["notes"] = pd.NA
    base = base[REFERENCE_COLUMNS]

    manual = manual.copy()
    manual["canonical_id"] = _normalize_series(manual["canonical_id"])
    manual["canonical_name"] = manual["canonical_name"].astype("string").str.strip()
    manual["iso3"] = _normalize_series(manual["iso3"])
    manual["iso2"] = manual["iso2"].astype("string").str.strip().str.upper()
    manual["entity_type"] = manual["entity_type"].astype("string").str.strip()
    manual["is_current"] = manual["is_current"].astype("string").str.strip()
    manual["valid_from"] = manual["valid_from"].astype("string").str.strip()
    manual["valid_to"] = manual["valid_to"].astype("string").str.strip()
    manual["world_bank_name"] = pd.NA
    manual["world_bank_region"] = pd.NA
    manual["world_bank_income_level"] = pd.NA
    manual["notes"] = manual["notes"].astype("string").str.strip()
    manual = manual[REFERENCE_COLUMNS]

    manual_ids = set(manual["canonical_id"].dropna())
    reference = pd.concat(
        [
            base.loc[~base["canonical_id"].isin(manual_ids)],
            manual,
        ],
        ignore_index=True,
    )
    reference = reference.sort_values(["canonical_id"]).reset_index(drop=True)

    missing_override_ids = sorted(set(overrides["canonical_id"].dropna()) - set(reference["canonical_id"].dropna()))
    if missing_override_ids:
        missing = ", ".join(missing_override_ids)
        raise ValueError(f"Unknown canonical ids in source overrides: {missing}")

    overrides = overrides.copy()
    overrides["source_id"] = overrides["source_id"].astype("string").str.strip()
    overrides["source_code"] = overrides["source_code"].astype("string").str.strip().str.upper()
    overrides["source_name"] = overrides["source_name"].astype("string").str.strip()
    overrides["canonical_id"] = _normalize_series(overrides["canonical_id"])
    overrides["match_priority"] = overrides["match_priority"].astype("string").str.strip()
    overrides["notes"] = overrides["notes"].astype("string").str.strip()
    overrides = overrides.merge(
        reference[
            [
                "canonical_id",
                "canonical_name",
                "iso3",
                "iso2",
                "entity_type",
                "is_current",
            ]
        ],
        on="canonical_id",
        how="left",
    )
    overrides = overrides[OVERRIDE_COLUMNS].sort_values(["source_id", "source_code", "source_name"]).reset_index(drop=True)

    reference_dir = settings.processed_dir / "reference"
    reference_dir.mkdir(parents=True, exist_ok=True)

    reference_path = reference_dir / "country_reference.csv"
    overrides_path = reference_dir / "source_country_overrides.csv"
    reference.to_csv(reference_path, index=False)
    overrides.to_csv(overrides_path, index=False)
    return [reference_path, overrides_path]


def load_country_reference(settings: Settings):
    pd = _require_pandas()
    path = settings.processed_dir / "reference" / "country_reference.csv"
    if not path.exists():
        raise FileNotFoundError(
            "Country reference is missing. Run `gravity-data build-country-reference` first."
        )
    return pd.read_csv(path, dtype="string")


def load_source_country_overrides(settings: Settings):
    pd = _require_pandas()
    path = settings.processed_dir / "reference" / "source_country_overrides.csv"
    if not path.exists():
        raise FileNotFoundError(
            "Source country overrides are missing. Run `gravity-data build-country-reference` first."
        )
    overrides = pd.read_csv(path, dtype="string")
    overrides["source_code_norm"] = _normalize_series(overrides["source_code"])
    overrides["source_name_norm"] = _normalize_series(overrides["source_name"])
    return overrides


def harmonize_source_frame(frame, settings: Settings, source_id: str, code_column: str | None = None, name_column: str | None = None):
    pd = _require_pandas()
    if code_column is None and name_column is None:
        raise ValueError("At least one of code_column or name_column must be provided.")

    reference = load_country_reference(settings).copy()
    overrides = load_source_country_overrides(settings)
    source_overrides = overrides.loc[overrides["source_id"] == source_id].copy()

    reference["canonical_id_norm"] = _normalize_series(reference["canonical_id"])
    reference["iso3_norm"] = _normalize_series(reference["iso3"])
    reference["canonical_name_norm"] = _normalize_series(reference["canonical_name"])

    working = frame.copy().reset_index(drop=True)
    working["canonical_id"] = pd.Series(pd.NA, index=working.index, dtype="string")
    working["match_method"] = pd.Series(pd.NA, index=working.index, dtype="string")

    if code_column is not None:
        working["source_code_norm"] = _normalize_series(working[code_column])

        code_override_map = source_overrides.loc[
            source_overrides["source_code_norm"].notna(),
            ["source_code_norm", "canonical_id"],
        ].drop_duplicates(subset=["source_code_norm"])
        if not code_override_map.empty:
            matched = working["source_code_norm"].map(code_override_map.set_index("source_code_norm")["canonical_id"])
            mask = working["canonical_id"].isna() & matched.notna()
            working.loc[mask, "canonical_id"] = matched[mask]
            working.loc[mask, "match_method"] = "override_code"

        direct_code_map = reference.loc[
            reference["canonical_id_norm"].notna(),
            ["canonical_id_norm", "canonical_id"],
        ].drop_duplicates(subset=["canonical_id_norm"])
        matched = working["source_code_norm"].map(direct_code_map.set_index("canonical_id_norm")["canonical_id"])
        mask = working["canonical_id"].isna() & matched.notna()
        working.loc[mask, "canonical_id"] = matched[mask]
        working.loc[mask, "match_method"] = "direct_canonical_code"

        iso3_map = reference.loc[
            reference["iso3_norm"].notna(),
            ["iso3_norm", "canonical_id"],
        ].drop_duplicates(subset=["iso3_norm"])
        matched = working["source_code_norm"].map(iso3_map.set_index("iso3_norm")["canonical_id"])
        mask = working["canonical_id"].isna() & matched.notna()
        working.loc[mask, "canonical_id"] = matched[mask]
        working.loc[mask, "match_method"] = "direct_iso3_code"

    if name_column is not None:
        working["source_name_norm"] = _normalize_series(working[name_column])

        name_override_map = source_overrides.loc[
            source_overrides["source_name_norm"].notna(),
            ["source_name_norm", "canonical_id"],
        ].drop_duplicates(subset=["source_name_norm"])
        if not name_override_map.empty:
            matched = working["source_name_norm"].map(name_override_map.set_index("source_name_norm")["canonical_id"])
            mask = working["canonical_id"].isna() & matched.notna()
            working.loc[mask, "canonical_id"] = matched[mask]
            working.loc[mask, "match_method"] = "override_name"

        direct_name_map = reference.loc[
            reference["canonical_name_norm"].notna(),
            ["canonical_name_norm", "canonical_id"],
        ].drop_duplicates(subset=["canonical_name_norm"])
        matched = working["source_name_norm"].map(direct_name_map.set_index("canonical_name_norm")["canonical_id"])
        mask = working["canonical_id"].isna() & matched.notna()
        working.loc[mask, "canonical_id"] = matched[mask]
        working.loc[mask, "match_method"] = "direct_name"

    working = working.merge(reference[REFERENCE_COLUMNS], on="canonical_id", how="left")
    unmatched = working.loc[working["canonical_id"].isna()].copy()
    return working, unmatched