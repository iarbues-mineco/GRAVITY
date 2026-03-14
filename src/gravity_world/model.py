from __future__ import annotations

import json
import math
import shutil
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

CEPII_STOCK_FEATURE_COLUMNS = [
    *CEPII_FEATURE_COLUMNS,
    "ln_migrant_stock_both_sexes_plus1",
]

CEPII_STOCK_UNEMPLOYMENT_FEATURE_COLUMNS = [
    *CEPII_STOCK_FEATURE_COLUMNS,
    "origin_unemployment_total_pct",
    "destination_unemployment_total_pct",
]

CEPII_STOCK_PPML_FEATURE_COLUMNS = [
    *CEPII_STOCK_FEATURE_COLUMNS,
]

DEFAULT_COMPARISON_COUNTRIES = ["ESP", "USA", "ITA", "FRA", "DEU"]
CURATED_OUTPUT_SUMMARY_NAME = "summary.md"


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
                "observed_balance_avg_annual": (observed_inflow_total_period - observed_outflow_total_period) / total_years,
                "fitted_balance_avg_annual_median": (
                    fitted_inflow_total_period_median - fitted_outflow_total_period_median
                ) / total_years,
                "fitted_balance_avg_annual_mean_smearing": (
                    fitted_inflow_total_period_mean_smearing - fitted_outflow_total_period_mean_smearing
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
                "observed_balance_total_period": observed_inflow_total_period - observed_outflow_total_period,
                "fitted_balance_total_period_median": fitted_inflow_total_period_median - fitted_outflow_total_period_median,
                "fitted_balance_total_period_mean_smearing": (
                    fitted_inflow_total_period_mean_smearing - fitted_outflow_total_period_mean_smearing
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
    spain_by_period,
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
            "Fitted values below use the median prediction `exp(Xb)`.",
            "",
            "| Country | Observed Inflow | Fitted Inflow Median | Observed Outflow | Fitted Outflow Median | Observed Balance | Fitted Balance Median |",
            "|---|---:|---:|---:|---:|---:|---:|",
        ]
    )

    for row in country_comparison.itertuples(index=False):
        lines.append(
            f"| {row.country_name} ({row.country_iso3}) | {row.observed_inflow_avg_annual:,.3f} | {row.fitted_inflow_avg_annual_median:,.3f} | {row.observed_outflow_avg_annual:,.3f} | {row.fitted_outflow_avg_annual_median:,.3f} | {row.observed_balance_avg_annual:,.3f} | {row.fitted_balance_avg_annual_median:,.3f} |"
        )

    lines.extend(
        [
            "",
            "## Spain By Period",
            "",
            "Observed and fitted values below are annualized flows. Fitted values use the median prediction `exp(Xb)`. Balance is `inflow - outflow`.",
            "",
            "| Period | Observed Inflow | Fitted Inflow Median | Observed Outflow | Fitted Outflow Median | Observed Balance | Fitted Balance Median |",
            "|---|---:|---:|---:|---:|---:|---:|",
        ]
    )

    for row in spain_by_period.itertuples(index=False):
        period_end_year = int(row.period_start_year) + int(row.period_length_years)
        period_label = f"{int(row.period_start_year)}-{period_end_year}"
        observed_balance = row.observed_inflow_spain - row.observed_outflow_spain
        fitted_balance_median = row.fitted_inflow_spain_median - row.fitted_outflow_spain_median
        lines.append(
            f"| {period_label} | {row.observed_inflow_spain:,.3f} | {row.fitted_inflow_spain_median:,.3f} | {row.observed_outflow_spain:,.3f} | {row.fitted_outflow_spain_median:,.3f} | {observed_balance:,.3f} | {fitted_balance_median:,.3f} |"
        )

    path.write_text("\n".join(lines), encoding="utf-8")

def _publish_curated_summary(settings: Settings, summary_md_path: Path) -> Path:
    output_dir = settings.root_dir / "output"
    output_dir.mkdir(parents=True, exist_ok=True)
    curated_summary_path = output_dir / CURATED_OUTPUT_SUMMARY_NAME
    shutil.copyfile(summary_md_path, curated_summary_path)
    return curated_summary_path


def _publish_named_summary(settings: Settings, summary_md_path: Path, output_name: str) -> Path:
    output_dir = settings.root_dir / "output"
    output_dir.mkdir(parents=True, exist_ok=True)
    curated_summary_path = output_dir / output_name
    shutil.copyfile(summary_md_path, curated_summary_path)
    return curated_summary_path


def _build_country_comparison_single(estimation, country_codes: list[str], fitted_total_period_column: str):
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
        fitted_inflow_total_period = float(estimation.loc[inflow_mask, fitted_total_period_column].sum())
        observed_outflow_total_period = float(estimation.loc[outflow_mask, "flow_total_period"].sum())
        fitted_outflow_total_period = float(estimation.loc[outflow_mask, fitted_total_period_column].sum())

        records.append(
            {
                "country_iso3": code,
                "country_name": country_names[code],
                "sample_years": total_years,
                "observed_inflow_avg_annual": observed_inflow_total_period / total_years,
                "fitted_inflow_avg_annual": fitted_inflow_total_period / total_years,
                "observed_outflow_avg_annual": observed_outflow_total_period / total_years,
                "fitted_outflow_avg_annual": fitted_outflow_total_period / total_years,
                "observed_balance_avg_annual": (observed_inflow_total_period - observed_outflow_total_period) / total_years,
                "fitted_balance_avg_annual": (fitted_inflow_total_period - fitted_outflow_total_period) / total_years,
                "observed_inflow_total_period": observed_inflow_total_period,
                "fitted_inflow_total_period": fitted_inflow_total_period,
                "observed_outflow_total_period": observed_outflow_total_period,
                "fitted_outflow_total_period": fitted_outflow_total_period,
                "observed_balance_total_period": observed_inflow_total_period - observed_outflow_total_period,
                "fitted_balance_total_period": fitted_inflow_total_period - fitted_outflow_total_period,
            }
        )

    return pd.DataFrame(records)


def _write_markdown_summary_ppml(
    path: Path,
    model_name: str,
    fit_stats: dict,
    coefficients,
    country_comparison,
    spain_by_period,
) -> None:
    lines = [
        f"# {model_name}",
        "",
        "## Fit Statistics",
        "",
        "| Metric | Value |",
        "|---|---:|",
        f"| Flow unit | {fit_stats['flow_unit']} |",
        f"| Prediction target | {fit_stats['level_prediction_mean']} |",
        f"| Observations in panel | {fit_stats['n_obs_total_panel']:,} |",
        f"| Observations used in estimation | {fit_stats['n_obs_estimation_sample']:,} |",
        f"| Zero-flow observations included | {fit_stats['n_obs_zero_flow_included']:,} |",
        f"| Parameters | {fit_stats['n_parameters']} |",
        f"| Converged | {fit_stats['converged']} |",
        f"| Iterations | {fit_stats['iterations']} |",
        f"| McFadden pseudo-R-squared | {fit_stats['pseudo_r_squared']:.6f} |",
        f"| RMSE (flow) | {fit_stats['rmse_flow']:.6f} |",
        f"| Deviance | {fit_stats['deviance']:.3f} |",
        f"| Poisson pseudo-log-likelihood | {fit_stats['log_likelihood']:.3f} |",
        f"| Null pseudo-log-likelihood | {fit_stats['null_log_likelihood']:.3f} |",
        f"| AIC | {fit_stats['aic']:.3f} |",
        f"| BIC | {fit_stats['bic']:.3f} |",
        f"| Mean observed flow | {fit_stats['mean_observed_flow']:.6f} |",
        f"| Mean fitted flow | {fit_stats['mean_fitted_flow']:.6f} |",
        "",
        "## Coefficients",
        "",
        "| Term | Estimate | Robust Std. Error | z-stat | p-value | 95% CI Low | 95% CI High |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ]

    for row in coefficients.itertuples(index=False):
        lines.append(
            f"| {row.term} | {row.estimate:.6f} | {row.std_error:.6f} | {row.z_stat:.6f} | {row.p_value_normal_approx:.6g} | {row.ci95_lower:.6f} | {row.ci95_upper:.6f} |"
        )

    lines.extend(
        [
            "",
            "## Average Annual Flow Comparison Over Full Sample",
            "",
            "Fitted values below are PPML fitted means.",
            "",
            "| Country | Observed Inflow | Fitted Inflow Mean | Observed Outflow | Fitted Outflow Mean | Observed Balance | Fitted Balance Mean |",
            "|---|---:|---:|---:|---:|---:|---:|",
        ]
    )

    for row in country_comparison.itertuples(index=False):
        lines.append(
            f"| {row.country_name} ({row.country_iso3}) | {row.observed_inflow_avg_annual:,.3f} | {row.fitted_inflow_avg_annual:,.3f} | {row.observed_outflow_avg_annual:,.3f} | {row.fitted_outflow_avg_annual:,.3f} | {row.observed_balance_avg_annual:,.3f} | {row.fitted_balance_avg_annual:,.3f} |"
        )

    lines.extend(
        [
            "",
            "## Spain By Period",
            "",
            "Observed and fitted values below are annualized flows. Fitted values are PPML fitted means. Balance is `inflow - outflow`.",
            "",
            "| Period | Observed Inflow | Fitted Inflow Mean | Observed Outflow | Fitted Outflow Mean | Observed Balance | Fitted Balance Mean |",
            "|---|---:|---:|---:|---:|---:|---:|",
        ]
    )

    for row in spain_by_period.itertuples(index=False):
        period_end_year = int(row.period_start_year) + int(row.period_length_years)
        period_label = f"{int(row.period_start_year)}-{period_end_year}"
        observed_balance = row.observed_inflow_spain - row.observed_outflow_spain
        fitted_balance = row.fitted_inflow_spain - row.fitted_outflow_spain
        lines.append(
            f"| {period_label} | {row.observed_inflow_spain:,.3f} | {row.fitted_inflow_spain:,.3f} | {row.observed_outflow_spain:,.3f} | {row.fitted_outflow_spain:,.3f} | {observed_balance:,.3f} | {fitted_balance:,.3f} |"
        )

    path.write_text("\n".join(lines), encoding="utf-8")


def _estimate_ppml_model(panel, feature_columns: list[str], model_name: str, output_prefix: str, output_dir: Path) -> list[Path]:
    pd = _require_pandas()
    import numpy as np

    required_columns = [
        'period_start_year',
        'period_end_year',
        'period_length_years',
        'origin_iso3',
        'origin_name',
        'destination_iso3',
        'destination_name',
        'flow',
        'flow_total_period',
        'flow_unit',
        *feature_columns,
    ]
    _validate_columns(panel, required_columns, f'Panel for {model_name}')

    estimation = panel.copy()
    estimation = estimation.replace([np.inf, -np.inf], np.nan)
    estimation = estimation.dropna(subset=['flow', *feature_columns]).reset_index(drop=True)
    estimation = estimation.loc[estimation['flow'].ge(0)].reset_index(drop=True)

    if estimation.empty:
        raise ValueError(f'No observations available for estimation in {model_name} after filtering.')

    X_df = _build_design_matrix(estimation, feature_columns)
    y = estimation['flow'].to_numpy(dtype=float)
    X = X_df.to_numpy(dtype=float)
    n_obs = int(len(y))
    n_params = int(X.shape[1])

    beta = np.zeros(n_params, dtype=float)
    mean_y = float(y.mean())
    beta[0] = math.log(mean_y if mean_y > 0 else 1e-12)

    converged = False
    iterations = 0
    max_iter = 100
    tolerance = 1e-8
    eta_clip = 30.0

    def loglike(current_beta):
        eta = np.clip(X @ current_beta, -eta_clip, eta_clip)
        mu = np.exp(eta)
        return float(np.sum(y * eta - mu))

    current_ll = loglike(beta)
    for iteration in range(1, max_iter + 1):
        eta = np.clip(X @ beta, -eta_clip, eta_clip)
        mu = np.exp(eta)
        xtwx = X.T @ (mu[:, None] * X)
        gradient = X.T @ (y - mu)
        step = np.linalg.pinv(xtwx) @ gradient

        step_scale = 1.0
        accepted = False
        while step_scale >= 1.0 / (2 ** 20):
            candidate_beta = beta + step_scale * step
            candidate_ll = loglike(candidate_beta)
            if candidate_ll >= current_ll - 1e-10:
                beta = candidate_beta
                current_ll = candidate_ll
                accepted = True
                break
            step_scale /= 2.0

        if not accepted:
            break

        iterations = iteration
        if float(np.max(np.abs(step_scale * step))) < tolerance:
            converged = True
            break

    eta = np.clip(X @ beta, -eta_clip, eta_clip)
    mu = np.exp(eta)
    xtwx = X.T @ (mu[:, None] * X)
    bread = np.linalg.pinv(xtwx)
    score = X * (y - mu)[:, None]
    meat = score.T @ score
    hc1_scale = n_obs / max(n_obs - n_params, 1)
    robust_cov = hc1_scale * (bread @ meat @ bread)
    std_errors = np.sqrt(np.clip(np.diag(robust_cov), a_min=0.0, a_max=None))
    z_stats = beta / std_errors
    p_values = np.array([_NormalApprox.two_sided_pvalue(float(value)) for value in z_stats])
    critical = _NormalApprox.critical_95()
    ci_lower = beta - critical * std_errors
    ci_upper = beta + critical * std_errors

    log_gamma = np.fromiter((math.lgamma(value + 1.0) for value in y), dtype=float, count=n_obs)
    null_mu = np.full(n_obs, y.mean(), dtype=float)
    null_ll = float(np.sum(y * np.log(np.clip(null_mu, 1e-12, None)) - null_mu - log_gamma))
    log_likelihood = float(np.sum(y * np.log(np.clip(mu, 1e-12, None)) - mu - log_gamma))
    pseudo_r_squared = float(1.0 - log_likelihood / null_ll) if null_ll != 0 else float('nan')
    deviance_terms = np.where(y > 0, y * np.log(np.clip(y / mu, 1e-12, None)) - (y - mu), - (y - mu))
    deviance = float(2.0 * np.sum(deviance_terms))
    aic = float(2 * n_params - 2 * log_likelihood)
    bic = float(np.log(n_obs) * n_params - 2 * log_likelihood)
    rmse_flow = float(np.sqrt(np.mean((y - mu) ** 2)))

    estimation['fitted_flow'] = mu
    estimation['raw_residual_flow'] = y - mu
    estimation['fitted_flow_total_period'] = mu * estimation['period_length_years'].astype(float)

    coefficients = pd.DataFrame(
        {
            'term': X_df.columns,
            'estimate': beta,
            'std_error': std_errors,
            'z_stat': z_stats,
            'p_value_normal_approx': p_values,
            'ci95_lower': ci_lower,
            'ci95_upper': ci_upper,
        }
    )

    fit_stats = {
        'model': model_name,
        'dependent_variable': 'flow',
        'flow_unit': str(estimation['flow_unit'].iloc[0]),
        'level_prediction_mean': 'E[flow|X] = exp(Xb)',
        'n_obs_total_panel': int(len(panel)),
        'n_obs_estimation_sample': n_obs,
        'n_obs_zero_flow_included': int((estimation['flow'] == 0).sum()),
        'n_parameters': n_params,
        'converged': converged,
        'iterations': iterations,
        'pseudo_r_squared': pseudo_r_squared,
        'rmse_flow': rmse_flow,
        'deviance': deviance,
        'log_likelihood': log_likelihood,
        'null_log_likelihood': null_ll,
        'aic': aic,
        'bic': bic,
        'mean_observed_flow': float(y.mean()),
        'mean_fitted_flow': float(mu.mean()),
        'feature_columns': feature_columns,
        'sample_years': float(estimation[['period_start_year', 'period_length_years']].drop_duplicates()['period_length_years'].sum()),
    }

    country_comparison = _build_country_comparison_single(estimation, DEFAULT_COMPARISON_COUNTRIES, 'fitted_flow_total_period')
    spain_mask_in = estimation['destination_iso3'].eq('ESP')
    spain_mask_out = estimation['origin_iso3'].eq('ESP')
    spain_by_period = (
        estimation.assign(
            observed_inflow_spain=np.where(spain_mask_in, estimation['flow'], 0.0),
            fitted_inflow_spain=np.where(spain_mask_in, estimation['fitted_flow'], 0.0),
            observed_outflow_spain=np.where(spain_mask_out, estimation['flow'], 0.0),
            fitted_outflow_spain=np.where(spain_mask_out, estimation['fitted_flow'], 0.0),
            observed_inflow_spain_total_period=np.where(spain_mask_in, estimation['flow_total_period'], 0.0),
            fitted_inflow_spain_total_period=np.where(spain_mask_in, estimation['fitted_flow_total_period'], 0.0),
            observed_outflow_spain_total_period=np.where(spain_mask_out, estimation['flow_total_period'], 0.0),
            fitted_outflow_spain_total_period=np.where(spain_mask_out, estimation['fitted_flow_total_period'], 0.0),
        )
        .groupby(['period_start_year', 'period_length_years'], as_index=False)[[
            'observed_inflow_spain',
            'fitted_inflow_spain',
            'observed_outflow_spain',
            'fitted_outflow_spain',
            'observed_inflow_spain_total_period',
            'fitted_inflow_spain_total_period',
            'observed_outflow_spain_total_period',
            'fitted_outflow_spain_total_period',
        ]]
        .sum()
    )
    spain_by_period['flow_unit'] = str(estimation['flow_unit'].iloc[0])

    coefficients_path = output_dir / f'{output_prefix}_coefficients.csv'
    fit_stats_path = output_dir / f'{output_prefix}_fit_stats.json'
    summary_path = output_dir / f'{output_prefix}_summary.txt'
    summary_md_path = output_dir / f'{output_prefix}_summary.md'
    fitted_sample_path = output_dir / f'{output_prefix}_fitted_sample.csv'
    spain_period_path = output_dir / f'{output_prefix}_spain_by_period.csv'
    country_comparison_path = output_dir / f'{output_prefix}_country_comparison.csv'

    coefficients.to_csv(coefficients_path, index=False)
    fit_stats_path.write_text(json.dumps(fit_stats, indent=2), encoding='utf-8')
    estimation.to_csv(fitted_sample_path, index=False)
    spain_by_period.to_csv(spain_period_path, index=False)
    country_comparison.to_csv(country_comparison_path, index=False)

    summary_lines = [
        f'{model_name} Summary',
        '',
        'Dependent variable: flow',
        f'Flow unit in panel and model outputs: {str(estimation["flow_unit"].iloc[0])}',
        'Prediction target: E[flow|X] = exp(Xb)',
        f'Observations in panel: {len(panel):,}',
        f'Observations used in estimation: {n_obs:,}',
        f'Zero-flow observations included: {int((estimation["flow"] == 0).sum()):,}',
        f'Parameters: {n_params}',
        f'Converged: {converged}',
        f'Iterations: {iterations}',
        f'McFadden pseudo-R-squared: {pseudo_r_squared:.6f}',
        f'RMSE (flow): {rmse_flow:.6f}',
        f'Deviance: {deviance:.3f}',
        f'Poisson pseudo-log-likelihood: {log_likelihood:.3f}',
        f'Null pseudo-log-likelihood: {null_ll:.3f}',
        f'AIC: {aic:.3f}',
        f'BIC: {bic:.3f}',
        f'Mean observed flow: {float(y.mean()):.6f}',
        f'Mean fitted flow: {float(mu.mean()):.6f}',
        '',
        'Coefficients:',
        'term,estimate,std_error,z_stat,p_value_normal_approx,ci95_lower,ci95_upper',
    ]
    for row in coefficients.itertuples(index=False):
        summary_lines.append(
            f'{row.term},{row.estimate:.6f},{row.std_error:.6f},{row.z_stat:.6f},{row.p_value_normal_approx:.6g},{row.ci95_lower:.6f},{row.ci95_upper:.6f}'
        )
    summary_path.write_text('\n'.join(summary_lines), encoding='utf-8')

    _write_markdown_summary_ppml(summary_md_path, model_name, fit_stats, coefficients, country_comparison, spain_by_period)

    return [
        summary_path,
        summary_md_path,
        coefficients_path,
        fit_stats_path,
        fitted_sample_path,
        spain_period_path,
        country_comparison_path,
    ]


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
        .groupby(["period_start_year", "period_length_years"], as_index=False)[[
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

    _write_markdown_summary(summary_md_path, model_name, fit_stats, coefficients, country_comparison, spain_totals, spain_by_period)

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


def estimate_cepii_stock_gravity_model(settings: Settings, output_dir: Path | None = None) -> list[Path]:
    pd = _require_pandas()
    panel_path = settings.processed_dir / "panels" / "bilateral_panel_cepii_stock.csv"
    if not panel_path.exists():
        raise FileNotFoundError(
            "CEPII stock panel is missing. Run `python -m gravity_world.cli assemble-cepii-stock-panel` first."
        )
    panel = pd.read_csv(panel_path)
    output_dir = output_dir or settings.processed_dir / "models"
    output_dir.mkdir(parents=True, exist_ok=True)
    output_paths = _estimate_model(
        panel,
        CEPII_STOCK_FEATURE_COLUMNS,
        "CEPII Gravity OLS With Migrant Stock",
        "cepii_stock_gravity_ols",
        output_dir,
    )
    summary_md_path = next(path for path in output_paths if path.name == "cepii_stock_gravity_ols_summary.md")
    curated_summary_path = _publish_curated_summary(settings, summary_md_path)
    return [*output_paths, curated_summary_path]


def estimate_cepii_stock_unemployment_gravity_model(settings: Settings, output_dir: Path | None = None) -> list[Path]:
    pd = _require_pandas()
    panel_path = settings.processed_dir / "panels" / "bilateral_panel_cepii_stock.csv"
    if not panel_path.exists():
        raise FileNotFoundError(
            "CEPII stock panel is missing. Run `python -m gravity_world.cli assemble-cepii-stock-panel` first."
        )
    panel = pd.read_csv(panel_path)
    output_dir = output_dir or settings.processed_dir / "models"
    output_dir.mkdir(parents=True, exist_ok=True)
    return _estimate_model(
        panel,
        CEPII_STOCK_UNEMPLOYMENT_FEATURE_COLUMNS,
        "CEPII Gravity OLS With Migrant Stock And Unemployment",
        "cepii_stock_unemployment_gravity_ols",
        output_dir,
    )


def estimate_cepii_stock_ppml_model(settings: Settings, output_dir: Path | None = None) -> list[Path]:
    pd = _require_pandas()
    panel_path = settings.processed_dir / "panels" / "bilateral_panel_cepii_stock.csv"
    if not panel_path.exists():
        raise FileNotFoundError(
            "CEPII stock panel is missing. Run `python -m gravity_world.cli assemble-cepii-stock-panel` first."
        )
    panel = pd.read_csv(panel_path)
    output_dir = output_dir or settings.processed_dir / "models"
    output_dir.mkdir(parents=True, exist_ok=True)
    output_paths = _estimate_ppml_model(
        panel,
        CEPII_STOCK_PPML_FEATURE_COLUMNS,
        "CEPII Gravity PPML With Migrant Stock",
        "cepii_stock_ppml",
        output_dir,
    )
    summary_md_path = next(path for path in output_paths if path.name == "cepii_stock_ppml_summary.md")
    curated_summary_path = _publish_named_summary(settings, summary_md_path, "ppml_summary.md")
    return [*output_paths, curated_summary_path]
