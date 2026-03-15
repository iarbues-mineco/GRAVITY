from __future__ import annotations

import json
import math
import shutil
from pathlib import Path

from .latent_geodesic import (
    EARTH_RADIUS_KM,
    _NormalApprox,
    _angles_to_unconstrained,
    _angles_to_vectors,
    _build_country_comparison,
    _build_reference_coordinates,
    _build_spain_importance_label_set,
    _compute_latent_distances,
    _compute_quadratic_displacement_penalty,
    _require_numpy,
    _require_pandas,
    _unconstrained_to_angles,
    _validate_columns,
    _write_map_svg,
)
from .settings import Settings

BASE_FEATURE_COLUMNS = [
    "ln_origin_population_total",
    "ln_destination_population_total",
    "ln_origin_gdp_pc_ppp_constant",
    "ln_destination_gdp_pc_ppp_constant",
    "contig",
    "comlang_off",
    "colony",
    "col45",
    "smctry",
    "ln_migrant_stock_both_sexes_plus1",
]


def _build_design_matrix_without_distance_fe(frame):
    pd = _require_pandas()
    design = frame[BASE_FEATURE_COLUMNS].copy()
    period_dummies = pd.get_dummies(frame["period_start_year"].astype("string"), prefix="period", dtype=float)
    if not period_dummies.empty:
        period_dummies = period_dummies.drop(columns=[sorted(period_dummies.columns)[0]])
    design = pd.concat([design, period_dummies], axis=1)
    return design, list(design.columns) + ["ln_latent_distance_km"]


def _build_initial_beta_fe(settings: Settings, term_names: list[str]):
    np = _require_numpy()
    beta = np.zeros(len(term_names), dtype=float)
    candidate_paths = [
        settings.processed_dir / "models" / "latent_geodesic_ppml_coefficients.csv",
        settings.processed_dir / "models" / "cepii_stock_ppml_coefficients.csv",
    ]
    pd = None
    for coefficients_path in candidate_paths:
        if not coefficients_path.exists():
            continue
        if pd is None:
            pd = _require_pandas()
        coefficients = pd.read_csv(coefficients_path)
        coef_map = dict(zip(coefficients["term"], coefficients["estimate"]))
        for idx, term in enumerate(term_names):
            source_term = term
            if term == "ln_latent_distance_km" and term not in coef_map:
                source_term = "ln_distance_km"
            if source_term in coef_map:
                beta[idx] = float(coef_map[source_term])
        break
    return beta


def _initialize_country_effects(y, origin_idx, destination_idx, n_countries: int):
    np = _require_numpy()
    origin_totals = np.bincount(origin_idx, weights=y, minlength=n_countries).astype(float)
    destination_totals = np.bincount(destination_idx, weights=y, minlength=n_countries).astype(float)
    mean_origin = max(float(origin_totals.mean()), 1e-12)
    mean_destination = max(float(destination_totals.mean()), 1e-12)
    origin_fe = np.log(np.clip(origin_totals / mean_origin, 1e-12, None))
    destination_fe = np.log(np.clip(destination_totals / mean_destination, 1e-12, None))
    shift = float(origin_fe.mean())
    origin_fe = origin_fe - shift
    destination_fe = destination_fe + shift
    return origin_fe, destination_fe, origin_totals, destination_totals


def _refit_country_effects(
    xbeta,
    y,
    origin_idx,
    destination_idx,
    origin_totals,
    destination_totals,
    origin_fe,
    destination_fe,
    max_iter: int = 80,
    tol: float = 1e-10,
):
    np = _require_numpy()
    eta_clip = 30.0
    origin_fe = origin_fe.copy()
    destination_fe = destination_fe.copy()
    iterations = 0
    converged = False

    for iteration in range(1, max_iter + 1):
        eta_origin = np.clip(xbeta + destination_fe[destination_idx], -eta_clip, eta_clip)
        base_origin = np.exp(eta_origin)
        denom_origin = np.bincount(origin_idx, weights=base_origin, minlength=len(origin_fe)).astype(float)
        delta_origin = np.zeros_like(origin_fe)
        mask_origin = (origin_totals > 0.0) & (denom_origin > 0.0)
        delta_origin[mask_origin] = np.log(origin_totals[mask_origin] / denom_origin[mask_origin])
        origin_fe = origin_fe + np.clip(delta_origin, -20.0, 20.0)

        eta_destination = np.clip(xbeta + origin_fe[origin_idx], -eta_clip, eta_clip)
        base_destination = np.exp(eta_destination)
        denom_destination = np.bincount(destination_idx, weights=base_destination, minlength=len(destination_fe)).astype(float)
        delta_destination = np.zeros_like(destination_fe)
        mask_destination = (destination_totals > 0.0) & (denom_destination > 0.0)
        delta_destination[mask_destination] = np.log(destination_totals[mask_destination] / denom_destination[mask_destination])
        destination_fe = destination_fe + np.clip(delta_destination, -20.0, 20.0)

        shift = float(origin_fe.mean())
        origin_fe = origin_fe - shift
        destination_fe = destination_fe + shift

        iterations = iteration
        max_update = max(
            float(np.max(np.abs(delta_origin[mask_origin]))) if np.any(mask_origin) else 0.0,
            float(np.max(np.abs(delta_destination[mask_destination]))) if np.any(mask_destination) else 0.0,
        )
        if max_update < tol:
            converged = True
            break

    eta = np.clip(xbeta + origin_fe[origin_idx] + destination_fe[destination_idx], -eta_clip, eta_clip)
    mu = np.exp(eta)
    objective = float((y * eta - mu).sum())
    return origin_fe, destination_fe, mu, objective, converged, iterations


def _fit_ppml_coefficients_with_country_fe(
    X,
    y,
    origin_idx,
    destination_idx,
    beta_init,
    origin_fe_init,
    destination_fe_init,
    origin_totals,
    destination_totals,
    max_iter: int = 60,
    tol: float = 1e-9,
):
    np = _require_numpy()
    beta = beta_init.astype(float).copy()
    origin_fe = origin_fe_init.astype(float).copy()
    destination_fe = destination_fe_init.astype(float).copy()

    xbeta = X @ beta
    origin_fe, destination_fe, mu, current_obj, fe_converged, fe_iterations = _refit_country_effects(
        xbeta,
        y,
        origin_idx,
        destination_idx,
        origin_totals,
        destination_totals,
        origin_fe,
        destination_fe,
    )

    converged = False
    iterations = 0
    fe_iterations_last = fe_iterations
    fe_converged_last = fe_converged

    for iteration in range(1, max_iter + 1):
        xtwx = X.T @ (mu[:, None] * X)
        gradient = X.T @ (y - mu)
        step = np.linalg.pinv(xtwx) @ gradient
        step_scale = 1.0
        accepted = False

        while step_scale >= 1.0 / (2**20):
            candidate_beta = beta + step_scale * step
            candidate_xbeta = X @ candidate_beta
            candidate_origin_fe, candidate_destination_fe, candidate_mu, candidate_obj, candidate_fe_converged, candidate_fe_iterations = _refit_country_effects(
                candidate_xbeta,
                y,
                origin_idx,
                destination_idx,
                origin_totals,
                destination_totals,
                origin_fe,
                destination_fe,
            )
            if candidate_obj >= current_obj - 1e-10:
                beta = candidate_beta
                origin_fe = candidate_origin_fe
                destination_fe = candidate_destination_fe
                mu = candidate_mu
                current_obj = candidate_obj
                fe_iterations_last = candidate_fe_iterations
                fe_converged_last = candidate_fe_converged
                accepted = True
                break
            step_scale /= 2.0

        if not accepted:
            break

        iterations = iteration
        if float(np.max(np.abs(step_scale * step))) < tol and fe_converged_last:
            converged = True
            break

    return beta, origin_fe, destination_fe, mu, current_obj, converged, iterations, fe_converged_last, fe_iterations_last


def _write_markdown_summary_fe(path: Path, fit_stats: dict, coefficients, country_comparison, spain_by_period) -> None:
    lines = [
        f"# {fit_stats['model_title']}",
        "",
        "## Fit Statistics",
        "",
        "| Metric | Value |",
        "|---|---:|",
        f"| Flow unit | {fit_stats['flow_unit']} |",
        f"| Prediction target | {fit_stats['prediction_target']} |",
        f"| Countries | {fit_stats['n_countries']} |",
        f"| Free country coordinates | {fit_stats['n_free_countries']} |",
        f"| Anchors | {fit_stats['anchors']} |",
        f"| Origin fixed effects | {fit_stats['origin_fe_count']} |",
        f"| Destination fixed effects | {fit_stats['destination_fe_count']} |",
        f"| Observations used in estimation | {fit_stats['n_obs_estimation_sample']:,} |",
        f"| Zero-flow observations included | {fit_stats['n_obs_zero_flow_included']:,} |",
        f"| Parameters | {fit_stats['n_parameters']} |",
        f"| Coordinate iterations | {fit_stats['coordinate_iterations']} |",
        f"| Slope iterations in last step | {fit_stats['beta_iterations_last_step']} |",
        f"| FE iterations in last step | {fit_stats['fe_iterations_last_step']} |",
        f"| Reference pseudo-log-likelihood | {fit_stats['reference_log_likelihood']:.3f} |",
        f"| Final pseudo-log-likelihood | {fit_stats['log_likelihood']:.3f} |",
        f"| Log-likelihood improvement | {fit_stats['log_likelihood_improvement']:.3f} |",
        f"| McFadden pseudo-R-squared | {fit_stats['pseudo_r_squared']:.6f} |",
        f"| RMSE (flow) | {fit_stats['rmse_flow']:.6f} |",
        f"| Mean observed flow | {fit_stats['mean_observed_flow']:.6f} |",
        f"| Mean fitted flow | {fit_stats['mean_fitted_flow']:.6f} |",
        f"| Mean displacement (km, free countries) | {fit_stats['mean_displacement_km']:.3f} |",
        f"| Median displacement (km, free countries) | {fit_stats['median_displacement_km']:.3f} |",
        f"| Max displacement (km, free countries) | {fit_stats['max_displacement_km']:.3f} |",
    ]
    if float(fit_stats.get("penalty_weight", 0.0)) > 0.0:
        lines.extend([
            f"| Quadratic penalty weight | {fit_stats['penalty_weight']:.6f} |",
            f"| Quadratic penalty value | {fit_stats['penalty_value']:.3f} |",
            f"| Penalized objective | {fit_stats['penalized_objective']:.3f} |",
            f"| Penalized objective improvement | {fit_stats['penalized_objective_improvement']:.3f} |",
        ])
    lines.extend([
        "",
        "## Slope Coefficients",
        "",
        "Standard errors below are conditional on the estimated country fixed effects.",
        "",
        "| Term | Estimate | Conditional Robust Std. Error | z-stat | p-value | 95% CI Low | 95% CI High |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ])
    for row in coefficients.itertuples(index=False):
        lines.append(
            f"| {row.term} | {row.estimate:.6f} | {row.std_error:.6f} | {row.z_stat:.6f} | {row.p_value_normal_approx:.6g} | {row.ci95_lower:.6f} | {row.ci95_upper:.6f} |"
        )
    lines.extend([
        "",
        "Full country fixed effects are written to the technical CSV outputs, not listed here.",
        "",
        "## Average Annual Flow Comparison Over Full Sample",
        "",
        "Fitted values below are latent-distance PPML fitted means.",
        "",
        "| Country | Observed Inflow | Fitted Inflow Mean | Observed Outflow | Fitted Outflow Mean | Observed Balance | Fitted Balance Mean |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ])
    for row in country_comparison.itertuples(index=False):
        lines.append(
            f"| {row.country_name} ({row.country_iso3}) | {row.observed_inflow_avg_annual:,.3f} | {row.fitted_inflow_avg_annual:,.3f} | {row.observed_outflow_avg_annual:,.3f} | {row.fitted_outflow_avg_annual:,.3f} | {row.observed_balance_avg_annual:,.3f} | {row.fitted_balance_avg_annual:,.3f} |"
        )
    lines.extend([
        "",
        "## Spain By Period",
        "",
        "Observed and fitted values below are annualized flows. Fitted values are latent-distance PPML fitted means. Balance is `inflow - outflow`.",
        "",
        "| Period | Observed Inflow | Fitted Inflow Mean | Observed Outflow | Fitted Outflow Mean | Observed Balance | Fitted Balance Mean |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ])
    for row in spain_by_period.itertuples(index=False):
        period_label = f"{int(row.period_start_year)}-{int(row.period_start_year) + int(row.period_length_years)}"
        observed_balance = row.observed_inflow_spain - row.observed_outflow_spain
        fitted_balance = row.fitted_inflow_spain - row.fitted_outflow_spain
        lines.append(
            f"| {period_label} | {row.observed_inflow_spain:,.3f} | {row.fitted_inflow_spain:,.3f} | {row.observed_outflow_spain:,.3f} | {row.fitted_outflow_spain:,.3f} | {observed_balance:,.3f} | {fitted_balance:,.3f} |"
        )
    path.write_text("\n".join(lines), encoding="utf-8")


def estimate_latent_geodesic_fe_ppml(
    settings: Settings,
    output_dir: Path | None = None,
    penalty_weight: float = 0.0,
    anchors: tuple[str, ...] | list[str] | None = None,
    progress_every: int = 1,
) -> list[Path]:
    pd = _require_pandas()
    np = _require_numpy()
    if penalty_weight < 0.0:
        raise ValueError("penalty_weight must be nonnegative.")

    anchor_codes = tuple(dict.fromkeys(code.upper() for code in (anchors or ("ESP", "FRA"))))
    if len(anchor_codes) < 2:
        raise ValueError("At least two anchor countries are required.")
    if progress_every < 1:
        raise ValueError("progress_every must be at least 1.")

    penalized = penalty_weight > 0.0
    model_prefix = "latent_geodesic_fe_penalized_ppml" if penalized else "latent_geodesic_fe_ppml"
    curated_base = "latent_geodesic_fe_penalized" if penalized else "latent_geodesic_fe"
    summary_title = (
        "Penalized Latent Geodesic PPML With Origin and Destination Fixed Effects"
        if penalized
        else "Latent Geodesic PPML With Origin and Destination Fixed Effects"
    )

    panel_path = settings.processed_dir / "panels" / "bilateral_panel_cepii_stock.csv"
    if not panel_path.exists():
        raise FileNotFoundError("CEPII stock panel is missing. Run `python -m gravity_world.cli assemble-cepii-stock-panel` first.")

    panel = pd.read_csv(panel_path)
    required_columns = [
        "period_start_year", "period_end_year", "period_length_years", "origin_canonical_id", "origin_iso3", "origin_name",
        "destination_canonical_id", "destination_iso3", "destination_name", "flow", "flow_total_period", "flow_unit", *BASE_FEATURE_COLUMNS,
    ]
    _validate_columns(panel, required_columns, "CEPII stock panel")

    estimation = panel.replace([np.inf, -np.inf], np.nan).dropna(subset=["flow", *BASE_FEATURE_COLUMNS]).copy()
    estimation = estimation.loc[estimation["flow"].ge(0)].reset_index(drop=True)
    if estimation.empty:
        raise ValueError("No observations available for latent geodesic FE PPML after filtering.")

    spain_label_codes = _build_spain_importance_label_set(panel)
    reference = _build_reference_coordinates(estimation, settings)
    missing_anchors = sorted(set(anchor_codes) - set(reference["iso3"]))
    if missing_anchors:
        raise ValueError("Anchor countries are not available in the estimation sample: " + ", ".join(missing_anchors))
    reference_by_canonical = reference.set_index("canonical_id")
    estimation["origin_country_index"] = estimation["origin_canonical_id"].map(reference_by_canonical["country_index"])
    estimation["destination_country_index"] = estimation["destination_canonical_id"].map(reference_by_canonical["country_index"])
    if estimation[["origin_country_index", "destination_country_index"]].isna().any().any():
        raise ValueError("Some panel countries could not be matched to the reference coordinate table.")

    design_base, term_names = _build_design_matrix_without_distance_fe(estimation)
    X_base = design_base.to_numpy(dtype=float)
    y = estimation["flow"].to_numpy(dtype=float)
    period_lengths = estimation["period_length_years"].to_numpy(dtype=float)
    origin_idx = estimation["origin_country_index"].to_numpy(dtype=int)
    destination_idx = estimation["destination_country_index"].to_numpy(dtype=int)

    reference_lat_rad = np.radians(reference["reference_lat_deg"].to_numpy(dtype=float))
    reference_lon_rad = np.radians(reference["reference_lon_deg"].to_numpy(dtype=float))
    country_count = len(reference)
    anchor_mask = reference["iso3"].isin(anchor_codes).to_numpy(dtype=bool)
    free_indices = np.flatnonzero(~anchor_mask)
    reference_free_vectors = _angles_to_vectors(reference_lat_rad[free_indices], reference_lon_rad[free_indices]) if len(free_indices) else None

    lat_param, lon_param = _angles_to_unconstrained(reference_lat_rad[free_indices], reference_lon_rad[free_indices])
    beta = _build_initial_beta_fe(settings, term_names)
    origin_fe, destination_fe, origin_totals, destination_totals = _initialize_country_effects(y, origin_idx, destination_idx, country_count)

    m_lat = np.zeros_like(lat_param)
    v_lat = np.zeros_like(lat_param)
    m_lon = np.zeros_like(lon_param)
    v_lon = np.zeros_like(lon_param)
    objective_history: list[float] = []
    best_state: dict[str, object] | None = None

    def compose_angles(current_lat_param, current_lon_param):
        lat = reference_lat_rad.copy()
        lon = reference_lon_rad.copy()
        lat_free, lon_free, dlat, dlon = _unconstrained_to_angles(current_lat_param, current_lon_param)
        lat[free_indices] = lat_free
        lon[free_indices] = lon_free
        return lat, lon, dlat, dlon

    print(
        "[latent-geodesic-fe] starting fit "
        f"| obs={len(y):,} countries={country_count} free={len(free_indices)} "
        f"anchors={','.join(anchor_codes)} penalty={penalty_weight:.4f}"
    )

    reference_log_distance, _, _, _, _, _ = _compute_latent_distances(reference_lat_rad, reference_lon_rad, origin_idx, destination_idx)
    reference_X = np.column_stack([X_base, reference_log_distance])
    _, _, _, _, reference_obj_no_const, _, _, _, _ = _fit_ppml_coefficients_with_country_fe(
        reference_X,
        y,
        origin_idx,
        destination_idx,
        beta.copy(),
        origin_fe.copy(),
        destination_fe.copy(),
        origin_totals,
        destination_totals,
    )

    print(f"[latent-geodesic-fe] reference fit complete | objective={reference_obj_no_const:.3f}")

    for iteration in range(1, 61):
        lat_rad, lon_rad, dlat_param, dlon_param = compose_angles(lat_param, lon_param)
        log_distance, distance_km, dlat_o, dlat_d, dlon_o, dlon_d = _compute_latent_distances(lat_rad, lon_rad, origin_idx, destination_idx)
        X = np.column_stack([X_base, log_distance])
        beta, origin_fe, destination_fe, mu, objective_no_const, beta_converged, beta_iterations, fe_converged, fe_iterations = _fit_ppml_coefficients_with_country_fe(
            X,
            y,
            origin_idx,
            destination_idx,
            beta,
            origin_fe,
            destination_fe,
            origin_totals,
            destination_totals,
        )
        penalty_value, penalty_grad_lat, penalty_grad_lon = _compute_quadratic_displacement_penalty(
            lat_rad, lon_rad, reference_lat_rad, reference_lon_rad, free_indices, penalty_weight, len(y)
        )
        penalized_objective = objective_no_const - penalty_value
        objective_history.append(penalized_objective)
        is_best = best_state is None or penalized_objective > float(best_state["objective"])
        if is_best:
            best_state = {
                "objective": float(penalized_objective),
                "objective_no_const": float(objective_no_const),
                "penalty_value": float(penalty_value),
                "beta": beta.copy(),
                "origin_fe": origin_fe.copy(),
                "destination_fe": destination_fe.copy(),
                "mu": mu.copy(),
                "log_distance": log_distance.copy(),
                "distance_km": distance_km.copy(),
                "lat_rad": lat_rad.copy(),
                "lon_rad": lon_rad.copy(),
                "iteration": iteration,
                "beta_converged": beta_converged,
                "beta_iterations": beta_iterations,
                "fe_converged": fe_converged,
                "fe_iterations": fe_iterations,
            }
        grad_common = (y - mu) * float(beta[-1])
        grad_lat = np.bincount(origin_idx, weights=grad_common * dlat_o, minlength=country_count) + np.bincount(destination_idx, weights=grad_common * dlat_d, minlength=country_count)
        grad_lon = np.bincount(origin_idx, weights=grad_common * dlon_o, minlength=country_count) + np.bincount(destination_idx, weights=grad_common * dlon_d, minlength=country_count)
        grad_lat_total = grad_lat + penalty_grad_lat
        grad_lon_total = grad_lon + penalty_grad_lon
        grad_lat_param = np.clip(grad_lat_total[free_indices] * dlat_param, -5000.0, 5000.0)
        grad_lon_param = np.clip(grad_lon_total[free_indices] * dlon_param, -5000.0, 5000.0)
        if len(free_indices):
            current_free_vectors = _angles_to_vectors(lat_rad[free_indices], lon_rad[free_indices])
            current_free_displacements = EARTH_RADIUS_KM * np.arccos(
                np.clip((reference_free_vectors * current_free_vectors).sum(axis=1), -1.0, 1.0)
            )
            mean_displacement_km = float(current_free_displacements.mean())
            max_displacement_km = float(current_free_displacements.max())
        else:
            mean_displacement_km = 0.0
            max_displacement_km = 0.0
        if iteration == 1 or iteration % progress_every == 0:
            delta_objective = float("nan") if len(objective_history) == 1 else penalized_objective - objective_history[-2]
            coord_grad = max(
                float(np.max(np.abs(grad_lat_param))) if len(grad_lat_param) else 0.0,
                float(np.max(np.abs(grad_lon_param))) if len(grad_lon_param) else 0.0,
            )
            print(
                "[latent-geodesic-fe] "
                f"iter={iteration:02d} objective={objective_no_const:.3f} penalized={penalized_objective:.3f} "
                f"delta={delta_objective:.3f} beta_iter={beta_iterations} fe_iter={fe_iterations} "
                f"mean_disp_km={mean_displacement_km:.1f} max_disp_km={max_displacement_km:.1f} "
                f"coord_grad={coord_grad:.3f} best={'yes' if is_best else 'no'}"
            )
        beta1 = 0.9
        beta2 = 0.999
        eps = 1e-8
        lr = 0.03
        m_lat = beta1 * m_lat + (1.0 - beta1) * grad_lat_param
        v_lat = beta2 * v_lat + (1.0 - beta2) * (grad_lat_param**2)
        m_lon = beta1 * m_lon + (1.0 - beta1) * grad_lon_param
        v_lon = beta2 * v_lon + (1.0 - beta2) * (grad_lon_param**2)
        lat_param = lat_param + lr * (m_lat / (1.0 - beta1**iteration)) / (np.sqrt(v_lat / (1.0 - beta2**iteration)) + eps)
        lon_param = lon_param + lr * (m_lon / (1.0 - beta1**iteration)) / (np.sqrt(v_lon / (1.0 - beta2**iteration)) + eps)
        if iteration >= 10 and len(objective_history) >= 6 and max(objective_history[-5:]) - min(objective_history[-5:]) < 1e-3:
            print("[latent-geodesic-fe] stopping: penalized objective stabilized over the last 5 iterations.")
            break

    if best_state is None:
        raise RuntimeError("Latent geodesic FE PPML failed to produce any evaluated state.")

    print(
        f"[latent-geodesic-fe] finished search | best_iter={int(best_state['iteration'])} "
        f"best_penalized={float(best_state['objective']):.3f} best_objective={float(best_state['objective_no_const']):.3f}"
    )

    final_beta = np.asarray(best_state["beta"], dtype=float)
    final_origin_fe = np.asarray(best_state["origin_fe"], dtype=float)
    final_destination_fe = np.asarray(best_state["destination_fe"], dtype=float)
    final_mu = np.asarray(best_state["mu"], dtype=float)
    final_log_distance = np.asarray(best_state["log_distance"], dtype=float)
    final_distance_km = np.asarray(best_state["distance_km"], dtype=float)
    final_lat_rad = np.asarray(best_state["lat_rad"], dtype=float)
    final_lon_rad = np.asarray(best_state["lon_rad"], dtype=float)
    final_iteration = int(best_state["iteration"])
    final_penalty_value = float(best_state["penalty_value"])
    final_X = np.column_stack([X_base, final_log_distance])

    score = final_X * (y - final_mu)[:, None]
    bread = np.linalg.pinv(final_X.T @ (final_mu[:, None] * final_X))
    meat = score.T @ score
    hc1_scale = len(y) / max(len(y) - len(final_beta), 1)
    robust_cov = hc1_scale * (bread @ meat @ bread)
    std_errors = np.sqrt(np.clip(np.diag(robust_cov), a_min=0.0, a_max=None))
    z_stats = np.divide(final_beta, std_errors, out=np.full_like(final_beta, np.nan), where=std_errors > 0)
    critical = _NormalApprox.critical_95()
    p_values = np.array([_NormalApprox.pvalue(float(value)) if math.isfinite(float(value)) else float("nan") for value in z_stats])
    ci_lower = final_beta - critical * std_errors
    ci_upper = final_beta + critical * std_errors

    log_gamma = np.fromiter((math.lgamma(value + 1.0) for value in y), dtype=float, count=len(y))
    null_mu = np.full(len(y), y.mean(), dtype=float)
    log_likelihood = float(np.sum(y * np.log(np.clip(final_mu, 1e-12, None)) - final_mu - log_gamma))
    null_log_likelihood = float(np.sum(y * np.log(np.clip(null_mu, 1e-12, None)) - null_mu - log_gamma))
    reference_log_likelihood = float(reference_obj_no_const - log_gamma.sum())
    penalized_objective = log_likelihood - final_penalty_value
    pseudo_r_squared = float(1.0 - log_likelihood / null_log_likelihood) if null_log_likelihood != 0 else float("nan")
    rmse_flow = float(np.sqrt(np.mean((y - final_mu) ** 2)))

    estimation = estimation.copy()
    estimation["latent_distance_km"] = final_distance_km
    estimation["ln_latent_distance_km"] = final_log_distance
    estimation["origin_fixed_effect"] = final_origin_fe[origin_idx]
    estimation["destination_fixed_effect"] = final_destination_fe[destination_idx]
    estimation["fitted_flow"] = final_mu
    estimation["fitted_flow_total_period"] = final_mu * period_lengths
    estimation["raw_residual_flow"] = y - final_mu

    slope_coefficients = pd.DataFrame({
        "term": term_names,
        "estimate": final_beta,
        "std_error": std_errors,
        "z_stat": z_stats,
        "p_value_normal_approx": p_values,
        "ci95_lower": ci_lower,
        "ci95_upper": ci_upper,
    })
    origin_effects = reference[["country_index", "iso3", "country_name"]].copy()
    origin_effects["origin_fixed_effect"] = final_origin_fe
    destination_effects = reference[["country_index", "iso3", "country_name"]].copy()
    destination_effects["destination_fixed_effect"] = final_destination_fe

    country_comparison = _build_country_comparison(estimation, "fitted_flow_total_period")
    spain_mask_in = estimation["destination_iso3"].eq("ESP")
    spain_mask_out = estimation["origin_iso3"].eq("ESP")
    spain_by_period = estimation.assign(
        observed_inflow_spain=np.where(spain_mask_in, estimation["flow"], 0.0),
        fitted_inflow_spain=np.where(spain_mask_in, estimation["fitted_flow"], 0.0),
        observed_outflow_spain=np.where(spain_mask_out, estimation["flow"], 0.0),
        fitted_outflow_spain=np.where(spain_mask_out, estimation["fitted_flow"], 0.0),
    ).groupby(["period_start_year", "period_length_years"], as_index=False)[["observed_inflow_spain", "fitted_inflow_spain", "observed_outflow_spain", "fitted_outflow_spain"]].sum()

    coordinates = reference.copy()
    coordinates["latent_lat_deg"] = np.degrees(final_lat_rad)
    coordinates["latent_lon_deg"] = np.degrees(final_lon_rad)
    ref_vec = _angles_to_vectors(np.radians(coordinates["reference_lat_deg"].to_numpy()), np.radians(coordinates["reference_lon_deg"].to_numpy()))
    lat_vec = _angles_to_vectors(np.radians(coordinates["latent_lat_deg"].to_numpy()), np.radians(coordinates["latent_lon_deg"].to_numpy()))
    coordinates["displacement_km"] = EARTH_RADIUS_KM * np.arccos(np.clip((ref_vec * lat_vec).sum(axis=1), -1.0, 1.0))
    coordinates["is_anchor"] = coordinates["iso3"].isin(anchor_codes)

    free_displacements = coordinates.loc[~coordinates["is_anchor"], "displacement_km"].to_numpy(dtype=float)
    fit_stats = {
        "model_title": summary_title,
        "flow_unit": str(estimation["flow_unit"].iloc[0]),
        "prediction_target": "E[flow|X] = exp(alpha_i + delta_j + Xb)",
        "n_countries": int(country_count),
        "n_free_countries": int(len(free_indices)),
        "anchors": ", ".join(anchor_codes),
        "origin_fe_count": int(country_count),
        "destination_fe_count": int(country_count),
        "n_obs_total_panel": int(len(panel)),
        "n_obs_estimation_sample": int(len(estimation)),
        "n_obs_zero_flow_included": int((estimation["flow"] == 0).sum()),
        "n_parameters": int(len(final_beta) + 2 * country_count + 2 * len(free_indices)),
        "coordinate_iterations": final_iteration,
        "beta_converged_last_step": bool(best_state["beta_converged"]),
        "beta_iterations_last_step": int(best_state["beta_iterations"]),
        "fe_converged_last_step": bool(best_state["fe_converged"]),
        "fe_iterations_last_step": int(best_state["fe_iterations"]),
        "reference_log_likelihood": reference_log_likelihood,
        "log_likelihood": log_likelihood,
        "log_likelihood_improvement": log_likelihood - reference_log_likelihood,
        "penalty_weight": float(penalty_weight),
        "penalty_value": final_penalty_value,
        "penalized_objective": penalized_objective,
        "penalized_objective_improvement": penalized_objective - reference_log_likelihood,
        "pseudo_r_squared": pseudo_r_squared,
        "rmse_flow": rmse_flow,
        "mean_observed_flow": float(y.mean()),
        "mean_fitted_flow": float(final_mu.mean()),
        "mean_displacement_km": float(free_displacements.mean()) if len(free_displacements) else 0.0,
        "median_displacement_km": float(np.median(free_displacements)) if len(free_displacements) else 0.0,
        "max_displacement_km": float(free_displacements.max()) if len(free_displacements) else 0.0,
    }

    output_dir = output_dir or settings.processed_dir / "models"
    output_dir.mkdir(parents=True, exist_ok=True)
    summary_txt_path = output_dir / f"{model_prefix}_summary.txt"
    summary_md_path = output_dir / f"{model_prefix}_summary.md"
    fit_stats_path = output_dir / f"{model_prefix}_fit_stats.json"
    coefficients_path = output_dir / f"{model_prefix}_coefficients.csv"
    origin_effects_path = output_dir / f"{model_prefix}_origin_effects.csv"
    destination_effects_path = output_dir / f"{model_prefix}_destination_effects.csv"
    fitted_sample_path = output_dir / f"{model_prefix}_fitted_sample.csv"
    spain_by_period_path = output_dir / f"{model_prefix}_spain_by_period.csv"
    country_comparison_path = output_dir / f"{model_prefix}_country_comparison.csv"
    coordinates_path = output_dir / f"{model_prefix}_country_coordinates.csv"
    map_path = output_dir / f"{model_prefix}_map.svg"

    summary_txt_lines = [
        f"{summary_title} Summary",
        "",
        f"Observations used in estimation: {len(estimation):,}",
        f"Countries: {country_count}",
        f"Free countries: {len(free_indices)}",
        f"Anchors: {fit_stats['anchors']}",
        f"Origin fixed effects: {country_count}",
        f"Destination fixed effects: {country_count}",
        f"Coordinate iterations: {final_iteration}",
        f"Final pseudo-log-likelihood: {log_likelihood:.3f}",
        f"Reference pseudo-log-likelihood: {reference_log_likelihood:.3f}",
        f"Log-likelihood improvement: {log_likelihood - reference_log_likelihood:.3f}",
        f"McFadden pseudo-R-squared: {pseudo_r_squared:.6f}",
        f"RMSE (flow): {rmse_flow:.6f}",
        f"Mean displacement (km, free countries): {fit_stats['mean_displacement_km']:.3f}",
        f"Median displacement (km, free countries): {fit_stats['median_displacement_km']:.3f}",
        f"Max displacement (km, free countries): {fit_stats['max_displacement_km']:.3f}",
    ]
    if penalized:
        summary_txt_lines.extend([
            f"Quadratic penalty weight: {penalty_weight:.6f}",
            f"Quadratic penalty value: {final_penalty_value:.3f}",
            f"Penalized objective: {penalized_objective:.3f}",
            f"Penalized objective improvement: {penalized_objective - reference_log_likelihood:.3f}",
        ])
    summary_txt_lines.extend([
        "",
        "Slope coefficients:",
        "term,estimate,std_error,z_stat,p_value_normal_approx,ci95_lower,ci95_upper",
        *[
            f"{row.term},{row.estimate:.6f},{row.std_error:.6f},{row.z_stat:.6f},{row.p_value_normal_approx:.6g},{row.ci95_lower:.6f},{row.ci95_upper:.6f}"
            for row in slope_coefficients.itertuples(index=False)
        ],
    ])

    summary_txt_path.write_text("\n".join(summary_txt_lines), encoding="utf-8")
    fit_stats_path.write_text(json.dumps(fit_stats, indent=2), encoding="utf-8")
    slope_coefficients.to_csv(coefficients_path, index=False)
    origin_effects.to_csv(origin_effects_path, index=False)
    destination_effects.to_csv(destination_effects_path, index=False)
    estimation.to_csv(fitted_sample_path, index=False)
    spain_by_period.to_csv(spain_by_period_path, index=False)
    country_comparison.to_csv(country_comparison_path, index=False)
    coordinates.to_csv(coordinates_path, index=False)
    _write_markdown_summary_fe(summary_md_path, fit_stats, slope_coefficients, country_comparison, spain_by_period)
    _write_map_svg(map_path, coordinates, anchor_codes=anchor_codes, spain_label_codes=spain_label_codes)

    curated_dir = settings.root_dir / "output"
    curated_dir.mkdir(parents=True, exist_ok=True)
    curated_summary_path = curated_dir / f"{curated_base}_summary.md"
    curated_map_path = curated_dir / f"{curated_base}_map.svg"
    shutil.copyfile(summary_md_path, curated_summary_path)
    shutil.copyfile(map_path, curated_map_path)

    return [
        summary_txt_path,
        summary_md_path,
        fit_stats_path,
        coefficients_path,
        origin_effects_path,
        destination_effects_path,
        fitted_sample_path,
        spain_by_period_path,
        country_comparison_path,
        coordinates_path,
        map_path,
        curated_summary_path,
        curated_map_path,
    ]
