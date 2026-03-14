from __future__ import annotations

import json
from pathlib import Path
from statistics import NormalDist

from .settings import Settings


REQUIRED_MODEL_COLUMNS = [
    "period_start_year",
    "origin_iso3",
    "destination_iso3",
    "flow",
    "flow_total_period",
    "flow_unit",
    "ln_flow",
    "ln_origin_population_total",
    "ln_destination_population_total",
    "ln_origin_gdp_pc_ppp_constant",
    "ln_destination_gdp_pc_ppp_constant",
]

BASE_FEATURE_COLUMNS = [
    "ln_origin_population_total",
    "ln_destination_population_total",
    "ln_origin_gdp_pc_ppp_constant",
    "ln_destination_gdp_pc_ppp_constant",
]


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


def _validate_columns(frame, required_columns: list[str]) -> None:
    missing = [column for column in required_columns if column not in frame.columns]
    if missing:
        actual = ", ".join(frame.columns)
        needed = ", ".join(missing)
        raise ValueError(f"Minimal panel is missing required columns: {needed}. Actual columns: {actual}")


def _build_design_matrix(frame):
    pd = _require_pandas()

    design = frame[BASE_FEATURE_COLUMNS].copy()
    period_dummies = pd.get_dummies(frame["period_start_year"].astype("string"), prefix="period", dtype=float)
    if not period_dummies.empty:
        baseline = sorted(period_dummies.columns)[0]
        period_dummies = period_dummies.drop(columns=[baseline])
    design = pd.concat([design, period_dummies], axis=1)
    design.insert(0, "const", 1.0)
    return design


def estimate_minimal_gravity_model(settings: Settings, output_dir: Path | None = None) -> list[Path]:
    pd = _require_pandas()
    import numpy as np

    panel_path = settings.processed_dir / "panels" / "bilateral_panel_minimal.csv"
    if not panel_path.exists():
        raise FileNotFoundError(
            "Minimal panel is missing. Run `python -m gravity_world.cli assemble-minimal-panel` first."
        )

    panel = pd.read_csv(panel_path)
    _validate_columns(panel, REQUIRED_MODEL_COLUMNS)

    estimation = panel.loc[panel["ln_flow"].notna()].copy()
    estimation = estimation.replace([np.inf, -np.inf], np.nan)
    estimation = estimation.dropna(subset=["ln_flow", *BASE_FEATURE_COLUMNS]).reset_index(drop=True)

    if estimation.empty:
        raise ValueError("No observations available for estimation after filtering to positive flows and complete regressors.")

    X_df = _build_design_matrix(estimation)
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
    fitted_flow = np.exp(fitted_ln) * smearing_factor
    fitted_flow_total_period = fitted_flow * estimation["period_length_years"].astype(float).to_numpy()

    estimation["fitted_ln_flow"] = fitted_ln
    estimation["residual_ln_flow"] = residuals
    estimation["fitted_flow"] = fitted_flow
    estimation["fitted_flow_total_period"] = fitted_flow_total_period

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
        "model": "minimal_gravity_ols_log_linear",
        "dependent_variable": "ln_flow",
        "flow_unit": str(estimation["flow_unit"].iloc[0]),
        "level_prediction": "exp(fitted_ln_flow) * smearing_factor",
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
    }

    spain_mask_in = estimation["destination_iso3"].eq("ESP")
    spain_mask_out = estimation["origin_iso3"].eq("ESP")
    spain_totals = pd.DataFrame(
        [
            {
                "scope": "all_periods_sum_of_annualized_flows",
                "flow_unit": str(estimation["flow_unit"].iloc[0]),
                "observed_inflow_spain": float(estimation.loc[spain_mask_in, "flow"].sum()),
                "fitted_inflow_spain": float(estimation.loc[spain_mask_in, "fitted_flow"].sum()),
                "observed_outflow_spain": float(estimation.loc[spain_mask_out, "flow"].sum()),
                "fitted_outflow_spain": float(estimation.loc[spain_mask_out, "fitted_flow"].sum()),
                "observed_total_spain": float(estimation.loc[spain_mask_in | spain_mask_out, "flow"].sum()),
                "fitted_total_spain": float(estimation.loc[spain_mask_in | spain_mask_out, "fitted_flow"].sum()),
                "observed_inflow_spain_total_period": float(estimation.loc[spain_mask_in, "flow_total_period"].sum()),
                "fitted_inflow_spain_total_period": float(estimation.loc[spain_mask_in, "fitted_flow_total_period"].sum()),
                "observed_outflow_spain_total_period": float(estimation.loc[spain_mask_out, "flow_total_period"].sum()),
                "fitted_outflow_spain_total_period": float(estimation.loc[spain_mask_out, "fitted_flow_total_period"].sum()),
            }
        ]
    )

    spain_by_period = (
        estimation.assign(
            observed_inflow_spain=np.where(spain_mask_in, estimation["flow"], 0.0),
            fitted_inflow_spain=np.where(spain_mask_in, estimation["fitted_flow"], 0.0),
            observed_outflow_spain=np.where(spain_mask_out, estimation["flow"], 0.0),
            fitted_outflow_spain=np.where(spain_mask_out, estimation["fitted_flow"], 0.0),
            observed_inflow_spain_total_period=np.where(spain_mask_in, estimation["flow_total_period"], 0.0),
            fitted_inflow_spain_total_period=np.where(spain_mask_in, estimation["fitted_flow_total_period"], 0.0),
            observed_outflow_spain_total_period=np.where(spain_mask_out, estimation["flow_total_period"], 0.0),
            fitted_outflow_spain_total_period=np.where(spain_mask_out, estimation["fitted_flow_total_period"], 0.0),
        )
        .groupby("period_start_year", as_index=False)[
            [
                "observed_inflow_spain",
                "fitted_inflow_spain",
                "observed_outflow_spain",
                "fitted_outflow_spain",
                "observed_inflow_spain_total_period",
                "fitted_inflow_spain_total_period",
                "observed_outflow_spain_total_period",
                "fitted_outflow_spain_total_period",
            ]
        ]
        .sum()
    )
    spain_by_period["flow_unit"] = str(estimation["flow_unit"].iloc[0])
    spain_by_period["observed_total_spain"] = spain_by_period["observed_inflow_spain"] + spain_by_period["observed_outflow_spain"]
    spain_by_period["fitted_total_spain"] = spain_by_period["fitted_inflow_spain"] + spain_by_period["fitted_outflow_spain"]
    spain_by_period["observed_total_spain_total_period"] = spain_by_period["observed_inflow_spain_total_period"] + spain_by_period["observed_outflow_spain_total_period"]
    spain_by_period["fitted_total_spain_total_period"] = spain_by_period["fitted_inflow_spain_total_period"] + spain_by_period["fitted_outflow_spain_total_period"]

    output_dir = output_dir or settings.processed_dir / "models"
    output_dir.mkdir(parents=True, exist_ok=True)

    coefficients_path = output_dir / "minimal_gravity_ols_coefficients.csv"
    fit_stats_path = output_dir / "minimal_gravity_ols_fit_stats.json"
    summary_path = output_dir / "minimal_gravity_ols_summary.txt"
    fitted_sample_path = output_dir / "minimal_gravity_ols_fitted_sample.csv"
    spain_totals_path = output_dir / "minimal_gravity_ols_spain_totals.csv"
    spain_period_path = output_dir / "minimal_gravity_ols_spain_by_period.csv"

    coefficients.to_csv(coefficients_path, index=False)
    fit_stats_path.write_text(json.dumps(fit_stats, indent=2), encoding="utf-8")
    estimation.to_csv(fitted_sample_path, index=False)
    spain_totals.to_csv(spain_totals_path, index=False)
    spain_by_period.to_csv(spain_period_path, index=False)

    summary_lines = [
        "Minimal Gravity OLS Summary",
        "",
        f"Dependent variable: ln_flow",
        f"Flow unit in panel and model outputs: {str(estimation['flow_unit'].iloc[0])}",
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
        "",
        "Coefficients:",
        "term,estimate,std_error,t_stat,p_value_normal_approx,ci95_lower,ci95_upper",
    ]
    for row in coefficients.itertuples(index=False):
        summary_lines.append(
            f"{row.term},{row.estimate:.6f},{row.std_error:.6f},{row.t_stat:.6f},{row.p_value_normal_approx:.6g},{row.ci95_lower:.6f},{row.ci95_upper:.6f}"
        )
    summary_lines.extend(
        [
            "",
            "Spain totals (estimation sample):",
            f"Observed inflow, annualized: {float(spain_totals.loc[0, 'observed_inflow_spain']):.6f}",
            f"Fitted inflow, annualized: {float(spain_totals.loc[0, 'fitted_inflow_spain']):.6f}",
            f"Observed outflow, annualized: {float(spain_totals.loc[0, 'observed_outflow_spain']):.6f}",
            f"Fitted outflow, annualized: {float(spain_totals.loc[0, 'fitted_outflow_spain']):.6f}",
            f"Observed total, annualized: {float(spain_totals.loc[0, 'observed_total_spain']):.6f}",
            f"Fitted total, annualized: {float(spain_totals.loc[0, 'fitted_total_spain']):.6f}",
            f"Observed inflow, period totals: {float(spain_totals.loc[0, 'observed_inflow_spain_total_period']):.6f}",
            f"Fitted inflow, period totals: {float(spain_totals.loc[0, 'fitted_inflow_spain_total_period']):.6f}",
            f"Observed outflow, period totals: {float(spain_totals.loc[0, 'observed_outflow_spain_total_period']):.6f}",
            f"Fitted outflow, period totals: {float(spain_totals.loc[0, 'fitted_outflow_spain_total_period']):.6f}",
        ]
    )
    summary_path.write_text("\n".join(summary_lines), encoding="utf-8")

    return [
        summary_path,
        coefficients_path,
        fit_stats_path,
        fitted_sample_path,
        spain_totals_path,
        spain_period_path,
    ]