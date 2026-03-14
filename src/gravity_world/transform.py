from __future__ import annotations

from pathlib import Path


def _require_pandas():
    try:
        import pandas as pd  # type: ignore
    except ImportError as exc:  # pragma: no cover
        raise RuntimeError("pandas is required for transformation commands.") from exc
    return pd


def assemble_country_year_covariates(raw_dir: Path, processed_dir: Path) -> Path:
    pd = _require_pandas()
    world_bank_dir = raw_dir / "world_bank"

    metadata_path = world_bank_dir / "country_metadata.csv"
    if not metadata_path.exists():
        raise FileNotFoundError(
            "World Bank country metadata is missing. Run `gravity-data download --source world_bank_wdi` first."
        )

    metadata = pd.read_csv(metadata_path, dtype={"id": "string", "iso2c": "string"})
    metadata = metadata.loc[metadata["region_id"].fillna("NA") != "NA", ["id", "name", "region_name", "income_level_name"]]
    metadata = metadata.rename(columns={"id": "country_iso2"})

    indicator_files = {
        "population_total": world_bank_dir / "population_total.csv",
        "gdp_pc_ppp_constant": world_bank_dir / "gdp_pc_ppp_constant.csv",
        "unemployment_total_pct": world_bank_dir / "unemployment_total_pct.csv",
        "gini_index": world_bank_dir / "gini_index.csv",
    }

    merged = None
    for slug, path in indicator_files.items():
        if not path.exists():
            raise FileNotFoundError(f"Missing indicator file: {path}")

        frame = pd.read_csv(path, dtype={"country_iso3": "string", "year": "string"})
        frame = frame[["country_iso3", "country_name", "year", "value"]].copy()
        frame["year"] = pd.to_numeric(frame["year"], errors="coerce").astype("Int64")
        frame = frame.rename(columns={"value": slug})

        if merged is None:
            merged = frame
        else:
            merged = merged.merge(frame[["country_iso3", "year", slug]], on=["country_iso3", "year"], how="outer")

    if merged is None:
        raise RuntimeError("No indicator files were loaded.")

    merged = merged.sort_values(["country_iso3", "year"]).reset_index(drop=True)
    output_path = processed_dir / "country_year_covariates.csv"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(output_path, index=False)
    return output_path
