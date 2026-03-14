from __future__ import annotations

import json
from pathlib import Path
from statistics import NormalDist

from .settings import Settings


COMMON_REQUIRED_COLUMNS = [
    "period_start_year",
    "period_length_years",
    "origin_iso3",
    "origin_name",
    "destination_iso3",
    "destination_name",
    "flow",
    "flow_total_period",
    "flow_unit",
    "ln_flow",
]

MINIMAL_FEATURE_COLUMNS = [
    "ln_origin_population_total",
    "ln_destination_population_total",
    "ln_origin_gdp_pc_ppp_constant",
    "ln_destination_gdp_pc_ppp_constant",
]

CEPII_FEATURE_COLUMNS = [
    *MINIMAL_FEATURE_COLUMNS,
    "ln_distance_km",
    "contig",
    "comlang_off",
    "colony",
    "col45",
    "smctry",
]

DEFAULT_COMPARISON_COUNTRIES = ["ESP", "USA", "ITA", "FRA", "DEU"]


class _NormalApprox:
    dist = NormalDist()

    @classmethod
    def two_sided_pvalue(cls, value: float) -> float:
        return 2.0 * (1.0 - cls.dist.cdf(abs(value)))

    @classmethod
    def critical_95(cls) -> float:
        return cls.dist.inv_cdf(0.975)


def _require_pandas():
    try:
        import pandas as pd  # type: ignore
    except ImportError as exc:  # pragma: no cover
        raise RuntimeError("pandas is required for model estimation commands.") from exc
    return pd


def _validate_columns(frame, required_columns: list[str], dataset_name: str) -> None:
    missing = [column for column in required_columns if column not in frame.columns]
    if missing:
        actual = ", ".join(frame.columns)
        needed = ", ".join(missing)
        raise ValueError(f"{dataset_name} is missing required columns: {needed}. Actual columns: {actual}")


def _build_design_matrix(frame, feature_columns: list[str]):
    pd = _require_pandas()

    design = frame[feature_columns].copy()
    period_dummies = pd.get_dummies(frame["period_start_year"].astype("string"), prefix="period", dtype=float)
    if not period_dummies.empty:
        baseline = sorted(period_dummies.columns)[0]
        period_dummies = period_dummies.drop(columns=[baseline])
    design = pd.concat([design, period_dummies], axis=1)
    design.insert(0, "const", 1.0)
    return design


def _extract_country_names(estimation, country_codes: list[str]) -> dict[str, str]:
    names: dict[str, str] = {}
    for code in country_codes:
        name_series = estimation.loc[estimation["origin_iso3"].eq(code), "origin_name"]
        if name_series.empty:
            name_series = estimation.loc[estimation["destination_iso3"].eq(code), "destination_name"]
        names[code] = str(name_series.iloc[0]) if not name_series.empty else code
    return names


def _build_country_comparison(estimation, country_codes: list[str]):
    pd = _require_pandas()

    unique_periods = (
        estimation[["period_start_year", "period_length_years"]]
        .drop_duplicates()
        .sort_values(["period_start_year"])
    )
    total_years = float(unique_periods["period_length_years"].sum())
    country_names = _extract_country_names(estimation, country_codes)

    records = []
    for code in country_codes:
        inflow_mask = estimation["destination_iso3"].eq(code)
        outflow_mask = estimation["origin_iso3"].eq(code)

        observed_inflow_total_period = float(estimation.loc[inflow_mask, "flow_total_period"].sum())
        fitted_inflow_total_period_median = float(estimation.loc[inflow_mask, "fitted_flow_total_period_median"].sum())
        fitted_inflow_total_period_mean_smearing = float(
            estimation.loc[inflow_mask, "fitted_flow_total_period_mean_smearing"].sum()
        )
        observed_outflow_total_period = float(estimation.loc[outflow_mask, "flow_total_period"].sum())
        fitted_outflow_total_period_median = float(estimation.loc[outflow_mask, "fitted_flow_total_period_median"].sum())
        fitted_outflow_total_period_mean_smearing = float(
            estimation.loc[outflow_mask, "fitted_flow_total_period_mean_smearing"].sum()
        )

        records.append(
            {
                "country_iso3": code,
                "country_name": country_names[code],
                "sample_years": total_years,
                "observed_inflow_avg_annual": observed_inflow_total_period / total_years,
                "fitted_inflow_avg_annual_median": fitted_inflow_total_period_median / total_years,
                "fitted_inflow_avg_annual_mean_smearing": fitted_inflow_total_period_mean_smearing / total_years,
                "observed_outflow_avg_annual": observed_outflow_total_period / total_years,
                "fitted_outflow_avg_annual_median": fitted_outflow_total_period_median / total_years,
                "fitted_outflow_avg_annual_mean_smearing": fitted_outflow_total_period_mean_smearing / total_years,
                "observed_total_avg_annual": (observed_inflow_total_period + observed_outflow_total_period) / total_years,
                "fitted_total_avg_annual_median": (
                    fitted_inflow_total_period_median + fitted_outflow_total_period_median
                ) / total_years,
                "fitted_total_avg_annual_mean_smearing": (
                    fitted_inflow_total_period_mean_smearing + fitted_outflow_total_period_mean_smearing
                ) / total_years,
                "observed_inflow_total_period": observed_inflow_total_period,
                "fitted_inflow_total_period_median": fitted_inflow_total_period_median,
                "fitted_inflow_total_period_mean_smearing": fitted_inflow_total_period_mean_smearing,
                "observed_outflow_total_period": observed_outflow_total_period,
                "fitted_outflow_total_period_median": fitted_outflow_total_period_median,
                "fitted_outflow_total_period_mean_smearing": fitted_outflow_total_period_mean_smearing,
                "observed_total_period": observed_inflow_total_period + observed_outflow_total_period,
                "fitted_total_period_median": fitted_inflow_total_period_median + fitted_outflow_total_period_median,
                "fitted_total_period_mean_smearing": (
                    fitted_inflow_total_period_mean_smearing + fitted_outflow_total_period_mean_smearing
                ),
            }
        )

    return pd.DataFrame(records)


def _write_markdown_summary(
    path: Path,
    model_name: str,
    fit_stats: dict,
    coefficients,
    country_comparison,
    spain_totals,
) -> None:
    lines = [
        f"# {model_name}",
        "",
        "## Fit Statistics",
        "",
        "| Metric | Value |",
        "|---|---:|",
        f"| Flow unit | {fit_stats['flow_unit']} |",
        f"| Level prediction (median) | `{fit_stats['median_level_prediction']}` |",
        f"| Level prediction (mean, smearing) | `{fit_stats['mean_level_prediction_smearing']}` |",
        f"| Observations in panel | {fit_stats['n_obs_total_panel']:,} |",
        f"| Observations used in estimation | {fit_stats['n_obs_estimation_sample']:,} |",
        f"| Zero-flow observations excluded | {fit_stats['n_obs_zero_flow_excluded']:,} |",
        f"| Parameters | {fit_stats['n_parameters']} |",
        f"| R-squared | {fit_stats['r_squared']:.6f} |",
        f"| Adjusted R-squared | {fit_stats['adj_r_squared']:.6f} |",
        f"| RMSE (ln flow) | {fit_stats['rmse_ln_flow']:.6f} |",
        f"| Residual std. error | {fit_stats['residual_std_error']:.6f} |",
        f"| Log-likelihood | {fit_stats['log_likelihood']:.3f} |",
        f"| AIC | {fit_stats['aic']:.3f} |",
        f"| BIC | {fit_stats['bic']:.3f} |",
        f"| Smearing factor | {fit_stats['smearing_factor']:.6f} |",
        "",
        "## Coefficients",
        "",
        "| Term | Estimate | Std. Error | t-stat | p-value | 95% CI Low | 95% CI High |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ]

    for row in coefficients.itertuples(index=False):
        lines.append(
            f"| {row.term} | {row.estimate:.6f} | {row.std_error:.6f} | {row.t_stat:.6f} | {row.p_value_normal_approx:.6g} | {row.ci95_lower:.6f} | {row.ci95_upper:.6f} |"
        )

    lines.extend(
        [
            "",
            "## Average Annual Flow Comparison Over Full Sample",
            "",
            "These averages are computed as total period flows over all periods divided by the total sample length in years.",
            "The two fitted columns are:",
            "- `Median`: `exp(Xb)`",
            "- `Mean (smearing)`: `exp(Xb) * mean(exp(residual))`",
            "",
            "| Country | Observed Inflow | Fitted Inflow Median | Fitted Inflow Mean | Observed Outflow | Fitted Outflow Median | Fitted Outflow Mean | Observed Total | Fitted Total Median | Fitted Total Mean |",
            "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
        ]
    )

    for row in country_comparison.itertuples(index=False):
        lines.append(
            f"| {row.country_name} ({row.country_iso3}) | {row.observed_inflow_avg_annual:,.3f} | {row.fitted_inflow_avg_annual_median:,.3f} | {row.fitted_inflow_avg_annual_mean_smearing:,.3f} | {row.observed_outflow_avg_annual:,.3f} | {row.fitted_outflow_avg_annual_median:,.3f} | {row.fitted_outflow_avg_annual_mean_smearing:,.3f} | {row.observed_total_avg_annual:,.3f} | {row.fitted_total_avg_annual_median:,.3f} | {row.fitted_total_avg_annual_mean_smearing:,.3f} |"
        )

    lines.extend(
        [
            "",
            "## Spain Totals",
            "",
            "| Measure | Value |",
            "|---|---:|",
            f"| Observed inflow, average annual | {float(spain_totals.loc[0, 'observed_inflow_avg_annual']):,.3f} |",
            f"| Fitted inflow, average annual, median | {float(spain_totals.loc[0, 'fitted_inflow_avg_annual_median']):,.3f} |",
            f"| Fitted inflow, average annual, mean (smearing) | {float(spain_totals.loc[0, 'fitted_inflow_avg_annual_mean_smearing']):,.3f} |",
            f"| Observed outflow, average annual | {float(spain_totals.loc[0, 'observed_outflow_avg_annual']):,.3f} |",
            f"| Fitted outflow, average annual, median | {float(spain_totals.loc[0, 'fitted_outflow_avg_annual_median']):,.3f} |",
            f"| Fitted outflow, average annual, mean (smearing) | {float(spain_totals.loc[0, 'fitted_outflow_avg_annual_mean_smearing']):,.3f} |",
            f"| Observed total, average annual | {float(spain_totals.loc[0, 'observed_total_avg_annual']):,.3f} |",
            f"| Fitted total, average annual, median | {float(spain_totals.loc[0, 'fitted_total_avg_annual_median']):,.3f} |",
            f"| Fitted total, average annual, mean (smearing) | {float(spain_totals.loc[0, 'fitted_total_avg_annual_mean_smearing']):,.3f} |",
            f"| Observed inflow, all period totals | {float(spain_totals.loc[0, 'observed_inflow_spain_total_period']):,.3f} |",
            f"| Fitted inflow, all period totals, median | {float(spain_totals.loc[0, 'fitted_inflow_spain_total_period_median']):,.3f} |",
            f"| Fitted inflow, all period totals, mean (smearing) | {float(spain_totals.loc[0, 'fitted_inflow_spain_total_period_mean_smearing']):,.3f} |",
            f"| Observed outflow, all period totals | {float(spain_totals.loc[0, 'observed_outflow_spain_total_period']):,.3f} |",
            f"| Fitted outflow, all period totals, median | {float(spain_totals.loc[0, 'fitted_outflow_spain_total_period_median']):,.3f} |",
            f"| Fitted outflow, all period totals, mean (smearing) | {float(spain_totals.loc[0, 'fitted_outflow_spain_total_period_mean_smearing']):,.3f} |",
        ]
    )

    path.write_text("\n".join(lines), encoding="utf-8")


def _estimate_model(panel, feature_columns: list[str], model_name: str, output_prefix: str, output_dir: Path) -> list[Path]:
    pd = _require_pandas()
    import numpy as np

    required_columns = [*COMMON_REQUIRED_COLUMNS, *feature_columns]
    _validate_columns(panel, required_columns, f"Panel for {model_name}")

    estimation = panel.loc[panel["ln_flow"].notna()].copy()
    estimation = estimation.replace([np.inf, -np.inf], np.nan)
    estimation = estimation.dropna(subset=["ln_flow", *feature_columns]).reset_index(drop=True)

    if estimation.empty:
        raise ValueError(f"No observations available for estimation in {model_name} after filtering.")

    X_df = _build_design_matrix(estimation, feature_columns)
    y = estimation["ln_flow"].to_numpy(dtype=float)
    X = X_df.to_numpy(dtype=float)

    beta, _, _, _ = np.linalg.lstsq(X, y, rcond=None)
    fitted_ln = X @ beta
    residuals = y - fitted_ln

    n_obs = int(len(y))
    n_params = int(X.shape[1])
    df_model = n_params - 1
    df_resid = n_obs - n_params
    ssr = float(np.dot(residuals, residuals))
    y_mean = float(y.mean())
    centered = y - y_mean
    sst = float(np.dot(centered, centered))
    r_squared = float(1.0 - ssr / sst) if sst > 0 else float("nan")
    adj_r_squared = float(1.0 - (1.0 - r_squared) * (n_obs - 1) / df_resid) if df_resid > 0 else float("nan")
    sigma2 = float(ssr / df_resid) if df_resid > 0 else float("nan")
    rmse = float((ssr / n_obs) ** 0.5)
    residual_std_error = float(sigma2 ** 0.5) if sigma2 >= 0 else float("nan")

    xtx_inv = np.linalg.pinv(X.T @ X)
    var_beta = sigma2 * xtx_inv
    std_errors = np.sqrt(np.clip(np.diag(var_beta), a_min=0.0, a_max=None))
    t_stats = beta / std_errors
    p_values = np.array([_NormalApprox.two_sided_pvalue(float(value)) for value in t_stats])
    critical = _NormalApprox.critical_95()
    ci_lower = beta - critical * std_errors
    ci_upper = beta + critical * std_errors

    sigma2_mle = float(ssr / n_obs)
    log_likelihood = float(-0.5 * n_obs * (np.log(2.0 * np.pi * sigma2_mle) + 1.0)) if sigma2_mle > 0 else float("nan")
    aic = float(2 * n_params - 2 * log_likelihood) if log_likelihood == log_likelihood else float("nan")
    bic = float(np.log(n_obs) * n_params - 2 * log_likelihood) if log_likelihood == log_likelihood else float("nan")

    smearing_factor = float(np.exp(residuals).mean())
    fitted_flow_median = np.exp(fitted_ln)
    fitted_flow_mean_smearing = fitted_flow_median * smearing_factor
    period_lengths = estimation["period_length_years"].astype(float).to_numpy()
    fitted_flow_total_period_median = fitted_flow_median * period_lengths
    fitted_flow_total_period_mean_smearing = fitted_flow_mean_smearing * period_lengths

    estimation["fitted_ln_flow"] = fitted_ln
    estimation["residual_ln_flow"] = residuals
    estimation["fitted_flow_median"] = fitted_flow_median
    estimation["fitted_flow_total_period_median"] = fitted_flow_total_period_median
    estimation["fitted_flow_mean_smearing"] = fitted_flow_mean_smearing
    estimation["fitted_flow_total_period_mean_smearing"] = fitted_flow_total_period_mean_smearing
    estimation["fitted_flow"] = fitted_flow_mean_smearing
    estimation["fitted_flow_total_period"] = fitted_flow_total_period_mean_smearing

    coefficients = pd.DataFrame(
        {
            "term": X_df.columns,
            "estimate": beta,
            "std_error": std_errors,
            "t_stat": t_stats,
            "p_value_normal_approx": p_values,
            "ci95_lower": ci_lower,
            "ci95_upper": ci_upper,
        }
    )

    fit_stats = {
        "model": model_name,
        "dependent_variable": "ln_flow",
        "flow_unit": str(estimation["flow_unit"].iloc[0]),
        "median_level_prediction": "exp(fitted_ln_flow)",
        "mean_level_prediction_smearing": "exp(fitted_ln_flow) * smearing_factor",
        "n_obs_total_panel": int(len(panel)),
        "n_obs_estimation_sample": n_obs,
        "n_obs_zero_flow_excluded": int((panel["flow"] <= 0).sum()),
        "n_parameters": n_params,
        "df_model": int(df_model),
        "df_resid": int(df_resid),
        "r_squared": r_squared,
        "adj_r_squared": adj_r_squared,
        "rmse_ln_flow": rmse,
        "residual_std_error": residual_std_error,
        "log_likelihood": log_likelihood,
        "aic": aic,
        "bic": bic,
        "smearing_factor": smearing_factor,
        "feature_columns": feature_columns,
        "sample_years": float(estimation[["period_start_year", "period_length_years"]].drop_duplicates()["period_length_years"].sum()),
    }

    country_comparison = _build_country_comparison(estimation, DEFAULT_COMPARISON_COUNTRIES)
    spain_totals = country_comparison.loc[country_comparison["country_iso3"].eq("ESP")].copy()
    spain_totals = spain_totals.rename(
        columns={
            "observed_inflow_avg_annual": "observed_inflow_avg_annual",
            "fitted_inflow_avg_annual_median": "fitted_inflow_avg_annual_median",
            "fitted_inflow_avg_annual_mean_smearing": "fitted_inflow_avg_annual_mean_smearing",
            "observed_outflow_avg_annual": "observed_outflow_avg_annual",
            "fitted_outflow_avg_annual_median": "fitted_outflow_avg_annual_median",
            "fitted_outflow_avg_annual_mean_smearing": "fitted_outflow_avg_annual_mean_smearing",
            "observed_total_avg_annual": "observed_total_avg_annual",
            "fitted_total_avg_annual_median": "fitted_total_avg_annual_median",
            "fitted_total_avg_annual_mean_smearing": "fitted_total_avg_annual_mean_smearing",
            "observed_inflow_total_period": "observed_inflow_spain_total_period",
            "fitted_inflow_total_period_median": "fitted_inflow_spain_total_period_median",
            "fitted_inflow_total_period_mean_smearing": "fitted_inflow_spain_total_period_mean_smearing",
            "observed_outflow_total_period": "observed_outflow_spain_total_period",
            "fitted_outflow_total_period_median": "fitted_outflow_spain_total_period_median",
            "fitted_outflow_total_period_mean_smearing": "fitted_outflow_spain_total_period_mean_smearing",
            "observed_total_period": "observed_total_spain_total_period",
            "fitted_total_period_median": "fitted_total_spain_total_period_median",
            "fitted_total_period_mean_smearing": "fitted_total_spain_total_period_mean_smearing",
        }
    )
    spain_totals.insert(0, "scope", "average_annual_over_full_sample")
    spain_totals.insert(1, "flow_unit", str(estimation["flow_unit"].iloc[0]))

    spain_mask_in = estimation["destination_iso3"].eq("ESP")
    spain_mask_out = estimation["origin_iso3"].eq("ESP")
    spain_by_period = (
        estimation.assign(
            observed_inflow_spain=np.where(spain_mask_in, estimation["flow"], 0.0),
            fitted_inflow_spain_median=np.where(spain_mask_in, estimation["fitted_flow_median"], 0.0),
            fitted_inflow_spain_mean_smearing=np.where(spain_mask_in, estimation["fitted_flow_mean_smearing"], 0.0),
            observed_outflow_spain=np.where(spain_mask_out, estimation["flow"], 0.0),
            fitted_outflow_spain_median=np.where(spain_mask_out, estimation["fitted_flow_median"], 0.0),
            fitted_outflow_spain_mean_smearing=np.where(spain_mask_out, estimation["fitted_flow_mean_smearing"], 0.0),
            observed_inflow_spain_total_period=np.where(spain_mask_in, estimation["flow_total_period"], 0.0),
            fitted_inflow_spain_total_period_median=np.where(
                spain_mask_in, estimation["fitted_flow_total_period_median"], 0.0
            ),
            fitted_inflow_spain_total_period_mean_smearing=np.where(
                spain_mask_in, estimation["fitted_flow_total_period_mean_smearing"], 0.0
            ),
            observed_outflow_spain_total_period=np.where(spain_mask_out, estimation["flow_total_period"], 0.0),
            fitted_outflow_spain_total_period_median=np.where(
                spain_mask_out, estimation["fitted_flow_total_period_median"], 0.0
            ),
            fitted_outflow_spain_total_period_mean_smearing=np.where(
                spain_mask_out, estimation["fitted_flow_total_period_mean_smearing"], 0.0
            ),
        )
        .groupby("period_start_year", as_index=False)[[
            "observed_inflow_spain",
            "fitted_inflow_spain_median",
            "fitted_inflow_spain_mean_smearing",
            "observed_outflow_spain",
            "fitted_outflow_spain_median",
            "fitted_outflow_spain_mean_smearing",
            "observed_inflow_spain_total_period",
            "fitted_inflow_spain_total_period_median",
            "fitted_inflow_spain_total_period_mean_smearing",
            "observed_outflow_spain_total_period",
            "fitted_outflow_spain_total_period_median",
            "fitted_outflow_spain_total_period_mean_smearing",
        ]]
        .sum()
    )
    spain_by_period["flow_unit"] = str(estimation["flow_unit"].iloc[0])
    spain_by_period["observed_total_spain"] = spain_by_period["observed_inflow_spain"] + spain_by_period["observed_outflow_spain"]
    spain_by_period["fitted_total_spain_median"] = (
        spain_by_period["fitted_inflow_spain_median"] + spain_by_period["fitted_outflow_spain_median"]
    )
    spain_by_period["fitted_total_spain_mean_smearing"] = (
        spain_by_period["fitted_inflow_spain_mean_smearing"] + spain_by_period["fitted_outflow_spain_mean_smearing"]
    )
    spain_by_period["observed_total_spain_total_period"] = (
        spain_by_period["observed_inflow_spain_total_period"] + spain_by_period["observed_outflow_spain_total_period"]
    )
    spain_by_period["fitted_total_spain_total_period_median"] = (
        spain_by_period["fitted_inflow_spain_total_period_median"]
        + spain_by_period["fitted_outflow_spain_total_period_median"]
    )
    spain_by_period["fitted_total_spain_total_period_mean_smearing"] = (
        spain_by_period["fitted_inflow_spain_total_period_mean_smearing"]
        + spain_by_period["fitted_outflow_spain_total_period_mean_smearing"]
    )

    coefficients_path = output_dir / f"{output_prefix}_coefficients.csv"
    fit_stats_path = output_dir / f"{output_prefix}_fit_stats.json"
    summary_path = output_dir / f"{output_prefix}_summary.txt"
    summary_md_path = output_dir / f"{output_prefix}_summary.md"
    fitted_sample_path = output_dir / f"{output_prefix}_fitted_sample.csv"
    spain_totals_path = output_dir / f"{output_prefix}_spain_totals.csv"
    spain_period_path = output_dir / f"{output_prefix}_spain_by_period.csv"
    country_comparison_path = output_dir / f"{output_prefix}_country_comparison.csv"

    coefficients.to_csv(coefficients_path, index=False)
    fit_stats_path.write_text(json.dumps(fit_stats, indent=2), encoding="utf-8")
    estimation.to_csv(fitted_sample_path, index=False)
    spain_totals.to_csv(spain_totals_path, index=False)
    spain_by_period.to_csv(spain_period_path, index=False)
    country_comparison.to_csv(country_comparison_path, index=False)

    summary_lines = [
        f"{model_name} Summary",
        "",
        "Dependent variable: ln_flow",
        f"Flow unit in panel and model outputs: {str(estimation['flow_unit'].iloc[0])}",
        "Level prediction (median): exp(fitted_ln_flow)",
        "Level prediction (mean, smearing): exp(fitted_ln_flow) * smearing_factor",
        f"Average annual comparison horizon: {fit_stats['sample_years']:.0f} years",
        f"Observations in panel: {len(panel):,}",
        f"Observations used in estimation: {n_obs:,}",
        f"Zero-flow observations excluded: {int((panel['flow'] <= 0).sum()):,}",
        f"Parameters: {n_params}",
        f"R-squared: {r_squared:.6f}",
        f"Adjusted R-squared: {adj_r_squared:.6f}",
        f"RMSE (ln flow): {rmse:.6f}",
        f"Residual std. error: {residual_std_error:.6f}",
        f"Log-likelihood: {log_likelihood:.3f}",
        f"AIC: {aic:.3f}",
        f"BIC: {bic:.3f}",
        f"Smearing factor: {smearing_factor:.6f}",
        f"Mean observed flow in estimation sample: {float(estimation['flow'].mean()):.6f}",
        f"Mean fitted flow, median prediction: {float(estimation['fitted_flow_median'].mean()):.6f}",
        f"Mean fitted flow, mean smearing prediction: {float(estimation['fitted_flow_mean_smearing'].mean()):.6f}",
        "",
        "Coefficients:",
        "term,estimate,std_error,t_stat,p_value_normal_approx,ci95_lower,ci95_upper",
    ]
    for row in coefficients.itertuples(index=False):
        summary_lines.append(
            f"{row.term},{row.estimate:.6f},{row.std_error:.6f},{row.t_stat:.6f},{row.p_value_normal_approx:.6g},{row.ci95_lower:.6f},{row.ci95_upper:.6f}"
        )
    summary_path.write_text("\n".join(summary_lines), encoding="utf-8")

    _write_markdown_summary(summary_md_path, model_name, fit_stats, coefficients, country_comparison, spain_totals)

    return [
        summary_path,
        summary_md_path,
        coefficients_path,
        fit_stats_path,
        fitted_sample_path,
        spain_totals_path,
        spain_period_path,
        country_comparison_path,
    ]


def estimate_minimal_gravity_model(settings: Settings, output_dir: Path | None = None) -> list[Path]:
    pd = _require_pandas()
    panel_path = settings.processed_dir / "panels" / "bilateral_panel_minimal.csv"
    if not panel_path.exists():
        raise FileNotFoundError(
            "Minimal panel is missing. Run `python -m gravity_world.cli assemble-minimal-panel` first."
        )
    panel = pd.read_csv(panel_path)
    output_dir = output_dir or settings.processed_dir / "models"
    output_dir.mkdir(parents=True, exist_ok=True)
    return _estimate_model(panel, MINIMAL_FEATURE_COLUMNS, "Minimal Gravity OLS", "minimal_gravity_ols", output_dir)


def estimate_cepii_gravity_model(settings: Settings, output_dir: Path | None = None) -> list[Path]:
    pd = _require_pandas()
    panel_path = settings.processed_dir / "panels" / "bilateral_panel_cepii.csv"
    if not panel_path.exists():
        raise FileNotFoundError(
            "CEPII panel is missing. Run `python -m gravity_world.cli assemble-cepii-panel` first."
        )
    panel = pd.read_csv(panel_path)
    output_dir = output_dir or settings.processed_dir / "models"
    output_dir.mkdir(parents=True, exist_ok=True)
    return _estimate_model(panel, CEPII_FEATURE_COLUMNS, "CEPII Gravity OLS", "cepii_gravity_ols", output_dir)
