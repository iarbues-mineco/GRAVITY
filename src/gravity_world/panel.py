from __future__ import annotations

from pathlib import Path

from .settings import Settings


REQUIRED_FLOW_COLUMNS = [
    "period_start_year",
    "period_end_year",
    "period_length_years",
    "period_label",
    "origin_canonical_id",
    "origin_iso3",
    "origin_name",
    "destination_canonical_id",
    "destination_iso3",
    "destination_name",
    "flow",
    "flow_total_period",
    "flow_measure",
    "flow_unit",
    "is_self_flow",
]

REQUIRED_COVARIATE_COLUMNS = [
    "country_iso3",
    "country_name",
    "year",
    "population_total",
    "gdp_pc_ppp_constant",
    "unemployment_total_pct",
    "gini_index",
]

REQUIRED_DYADIC_COLUMNS = [
    "origin_canonical_id",
    "destination_canonical_id",
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
]

REQUIRED_STOCK_COLUMNS = [
    "stock_year",
    "origin_canonical_id",
    "destination_canonical_id",
    "migrant_stock_both_sexes",
    "migrant_stock_male",
    "migrant_stock_female",
]


def _require_pandas():
    try:
        import pandas as pd  # type: ignore
    except ImportError as exc:  # pragma: no cover
        raise RuntimeError("pandas is required for panel assembly commands.") from exc
    return pd


def _validate_columns(frame, required_columns: list[str], dataset_name: str) -> None:
    missing = [column for column in required_columns if column not in frame.columns]
    if missing:
        actual = ", ".join(frame.columns)
        needed = ", ".join(missing)
        raise ValueError(f"{dataset_name} is missing required columns: {needed}. Actual columns: {actual}")


def _prepare_covariates(frame):
    covariates = frame.copy()
    covariates = covariates.loc[covariates["country_iso3"].notna() & covariates["year"].notna()].copy()
    covariates = covariates.drop_duplicates(subset=["country_iso3", "year"], keep="first")
    return covariates


def _prefix_covariates(frame, prefix: str):
    renamed = frame.rename(
        columns={
            "country_iso3": f"{prefix}_canonical_id",
            "country_name": f"{prefix}_covariate_name",
            "year": "period_start_year",
            "population_total": f"{prefix}_population_total",
            "gdp_pc_ppp_constant": f"{prefix}_gdp_pc_ppp_constant",
            "unemployment_total_pct": f"{prefix}_unemployment_total_pct",
            "gini_index": f"{prefix}_gini_index",
        }
    )
    return renamed


def _derive_minimal_drop_reason(frame):
    pd = _require_pandas()

    def build_reason(row):
        reasons = []
        if bool(row.get("is_self_flow")):
            reasons.append("self_flow")
        if pd.isna(row.get("origin_canonical_id")):
            reasons.append("unmatched_origin")
        if pd.isna(row.get("destination_canonical_id")):
            reasons.append("unmatched_destination")
        if pd.isna(row.get("origin_population_total")) or pd.isna(row.get("origin_gdp_pc_ppp_constant")):
            reasons.append("missing_origin_minimal_covariates")
        if pd.isna(row.get("destination_population_total")) or pd.isna(row.get("destination_gdp_pc_ppp_constant")):
            reasons.append("missing_destination_minimal_covariates")
        return ";".join(reasons)

    return frame.apply(build_reason, axis=1)


def _add_log_columns(frame):
    import numpy as np

    for column in [
        "flow",
        "origin_population_total",
        "destination_population_total",
        "origin_gdp_pc_ppp_constant",
        "destination_gdp_pc_ppp_constant",
        "distance_km",
    ]:
        if column not in frame.columns:
            continue
        values = frame[column].astype(float)
        logged = values.copy()
        logged[:] = np.nan
        positive_mask = values > 0
        logged.loc[positive_mask] = np.log(values.loc[positive_mask])
        frame[f"ln_{column}"] = logged

    if "migrant_stock_both_sexes" in frame.columns:
        stock = frame["migrant_stock_both_sexes"].astype(float)
        frame["ln_migrant_stock_both_sexes_plus1"] = np.log1p(stock.clip(lower=0))
    return frame


def assemble_minimal_bilateral_panel(settings: Settings, output_dir: Path | None = None) -> list[Path]:
    pd = _require_pandas()

    flows_path = settings.processed_dir / "flows" / "bilateral_flows_5y_abel_cohen.csv"
    covariates_path = settings.processed_dir / "country_year_covariates.csv"

    if not flows_path.exists():
        raise FileNotFoundError(
            "Normalized Abel flow file is missing. Run `python -m gravity_world.cli normalize-abel-cohen` first."
        )
    if not covariates_path.exists():
        raise FileNotFoundError(
            "Country-year covariates are missing. Run `python -m gravity_world.cli assemble-country-year` first."
        )

    flows = pd.read_csv(flows_path)
    covariates = pd.read_csv(covariates_path)
    _validate_columns(flows, REQUIRED_FLOW_COLUMNS, "Normalized flow file")
    _validate_columns(covariates, REQUIRED_COVARIATE_COLUMNS, "Country-year covariates")

    covariates = _prepare_covariates(covariates)
    origin_covariates = _prefix_covariates(covariates, "origin")
    destination_covariates = _prefix_covariates(covariates, "destination")

    merged = flows.merge(origin_covariates, on=["origin_canonical_id", "period_start_year"], how="left")
    merged = merged.merge(destination_covariates, on=["destination_canonical_id", "period_start_year"], how="left")
    merged["covariate_year"] = merged["period_start_year"]
    merged["drop_reason"] = _derive_minimal_drop_reason(merged)

    minimal_panel = merged.loc[merged["drop_reason"].eq("")].copy()
    minimal_panel = _add_log_columns(minimal_panel)

    output_dir = output_dir or settings.processed_dir / "panels"
    output_dir.mkdir(parents=True, exist_ok=True)

    panel_path = output_dir / "bilateral_panel_minimal.csv"
    dropped_path = output_dir / "bilateral_panel_minimal_dropped_rows.csv"

    panel_columns = [
        "source_dataset",
        "period_start_year",
        "period_end_year",
        "period_length_years",
        "period_label",
        "covariate_year",
        "origin_canonical_id",
        "origin_iso3",
        "origin_name",
        "destination_canonical_id",
        "destination_iso3",
        "destination_name",
        "flow",
        "flow_total_period",
        "flow_measure",
        "flow_unit",
        "origin_population_total",
        "destination_population_total",
        "origin_gdp_pc_ppp_constant",
        "destination_gdp_pc_ppp_constant",
        "origin_unemployment_total_pct",
        "destination_unemployment_total_pct",
        "origin_gini_index",
        "destination_gini_index",
        "ln_flow",
        "ln_origin_population_total",
        "ln_destination_population_total",
        "ln_origin_gdp_pc_ppp_constant",
        "ln_destination_gdp_pc_ppp_constant",
    ]
    minimal_panel = minimal_panel[panel_columns].sort_values(
        ["period_start_year", "origin_canonical_id", "destination_canonical_id"]
    ).reset_index(drop=True)

    dropped_columns = [
        "source_dataset",
        "period_start_year",
        "period_end_year",
        "period_length_years",
        "period_label",
        "origin_code_raw",
        "origin_canonical_id",
        "destination_code_raw",
        "destination_canonical_id",
        "flow",
        "flow_total_period",
        "drop_reason",
    ]
    dropped_rows = merged.loc[merged["drop_reason"].ne(""), dropped_columns].copy()
    dropped_rows = dropped_rows.sort_values(
        ["period_start_year", "origin_code_raw", "destination_code_raw"]
    ).reset_index(drop=True)

    minimal_panel.to_csv(panel_path, index=False)
    dropped_rows.to_csv(dropped_path, index=False)
    return [panel_path, dropped_path]


def assemble_cepii_bilateral_panel(settings: Settings, output_dir: Path | None = None) -> list[Path]:
    pd = _require_pandas()

    minimal_panel_path = settings.processed_dir / "panels" / "bilateral_panel_minimal.csv"
    dyadic_path = settings.processed_dir / "dyadic" / "cepii_geodist_controls.csv"

    if not minimal_panel_path.exists():
        raise FileNotFoundError(
            "Minimal panel is missing. Run `python -m gravity_world.cli assemble-minimal-panel` first."
        )
    if not dyadic_path.exists():
        raise FileNotFoundError(
            "CEPII dyadic controls are missing. Run `python -m gravity_world.cli normalize-cepii` first."
        )

    minimal_panel = pd.read_csv(minimal_panel_path)
    dyadic = pd.read_csv(dyadic_path)
    _validate_columns(dyadic, REQUIRED_DYADIC_COLUMNS, "CEPII dyadic controls")

    dyadic = dyadic[
        [
            "origin_canonical_id",
            "destination_canonical_id",
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
        ]
    ].drop_duplicates(subset=["origin_canonical_id", "destination_canonical_id"], keep="first")
    merged = minimal_panel.merge(dyadic, on=["origin_canonical_id", "destination_canonical_id"], how="left")

    merged["cepii_drop_reason"] = ""
    missing_distance_mask = merged["distance_km"].isna()
    merged.loc[missing_distance_mask, "cepii_drop_reason"] = "missing_cepii_distance"

    extended_panel = merged.loc[merged["cepii_drop_reason"].eq("")].copy()
    extended_panel = _add_log_columns(extended_panel)

    output_dir = output_dir or settings.processed_dir / "panels"
    output_dir.mkdir(parents=True, exist_ok=True)

    panel_path = output_dir / "bilateral_panel_cepii.csv"
    dropped_path = output_dir / "bilateral_panel_cepii_dropped_rows.csv"

    panel_columns = list(minimal_panel.columns) + [
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
    ]
    extended_panel = extended_panel[panel_columns].sort_values(
        ["period_start_year", "origin_canonical_id", "destination_canonical_id"]
    ).reset_index(drop=True)

    dropped_columns = [
        "period_start_year",
        "origin_canonical_id",
        "destination_canonical_id",
        "flow",
        "flow_total_period",
        "cepii_drop_reason",
    ]
    dropped_rows = merged.loc[merged["cepii_drop_reason"].ne(""), dropped_columns].copy()
    dropped_rows = dropped_rows.sort_values(
        ["period_start_year", "origin_canonical_id", "destination_canonical_id"]
    ).reset_index(drop=True)

    extended_panel.to_csv(panel_path, index=False)
    dropped_rows.to_csv(dropped_path, index=False)
    return [panel_path, dropped_path]


def assemble_cepii_stock_bilateral_panel(settings: Settings, output_dir: Path | None = None) -> list[Path]:
    pd = _require_pandas()

    cepii_panel_path = settings.processed_dir / "panels" / "bilateral_panel_cepii.csv"
    stock_path = settings.processed_dir / "stocks" / "bilateral_migrant_stock_un_desa.csv"

    if not cepii_panel_path.exists():
        raise FileNotFoundError(
            "CEPII panel is missing. Run `python -m gravity_world.cli assemble-cepii-panel` first."
        )
    if not stock_path.exists():
        raise FileNotFoundError(
            "UN DESA stock file is missing. Run `python -m gravity_world.cli normalize-un-desa-stock` first."
        )

    cepii_panel = pd.read_csv(cepii_panel_path)
    stock = pd.read_csv(stock_path)
    _validate_columns(stock, REQUIRED_STOCK_COLUMNS, "UN DESA bilateral stock file")

    stock = stock.copy()
    stock["stock_year"] = pd.to_numeric(stock["stock_year"], errors="coerce").astype("Int64")
    stock = stock.rename(columns={"stock_year": "period_start_year"})
    stock = stock.drop_duplicates(
        subset=["period_start_year", "origin_canonical_id", "destination_canonical_id"],
        keep="first",
    )

    stock_columns = [
        "period_start_year",
        "origin_canonical_id",
        "destination_canonical_id",
        "stock_reference",
        "stock_unit",
        "migrant_stock_both_sexes",
        "migrant_stock_male",
        "migrant_stock_female",
    ]
    merged = cepii_panel.merge(
        stock[stock_columns],
        on=["period_start_year", "origin_canonical_id", "destination_canonical_id"],
        how="left",
    )

    merged["stock_drop_reason"] = ""
    missing_stock_mask = merged["migrant_stock_both_sexes"].isna()
    merged.loc[missing_stock_mask, "stock_drop_reason"] = "missing_un_desa_stock"

    stock_panel = merged.loc[merged["stock_drop_reason"].eq("")].copy()
    stock_panel = _add_log_columns(stock_panel)

    output_dir = output_dir or settings.processed_dir / "panels"
    output_dir.mkdir(parents=True, exist_ok=True)

    panel_path = output_dir / "bilateral_panel_cepii_stock.csv"
    dropped_path = output_dir / "bilateral_panel_cepii_stock_dropped_rows.csv"

    panel_columns = list(cepii_panel.columns) + [
        "stock_reference",
        "stock_unit",
        "migrant_stock_both_sexes",
        "migrant_stock_male",
        "migrant_stock_female",
        "ln_migrant_stock_both_sexes_plus1",
    ]
    stock_panel = stock_panel[panel_columns].sort_values(
        ["period_start_year", "origin_canonical_id", "destination_canonical_id"]
    ).reset_index(drop=True)

    dropped_columns = [
        "period_start_year",
        "origin_canonical_id",
        "destination_canonical_id",
        "flow",
        "flow_total_period",
        "stock_drop_reason",
    ]
    dropped_rows = merged.loc[merged["stock_drop_reason"].ne(""), dropped_columns].copy()
    dropped_rows = dropped_rows.sort_values(
        ["period_start_year", "origin_canonical_id", "destination_canonical_id"]
    ).reset_index(drop=True)

    stock_panel.to_csv(panel_path, index=False)
    dropped_rows.to_csv(dropped_path, index=False)
    return [panel_path, dropped_path]
