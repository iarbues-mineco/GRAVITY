from __future__ import annotations

import json
import math
import shutil
from pathlib import Path

from .settings import Settings

EARTH_RADIUS_KM = 6371.0088
MIN_DISTANCE_KM = 50.0
HALF_PI = math.pi / 2.0
DEFAULT_COMPARISON_COUNTRIES = ["ESP", "USA", "ITA", "FRA", "DEU"]
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
PENALTY_DISTANCE_SCALE_KM = 1000.0


class _NormalApprox:
    @staticmethod
    def cdf(value: float) -> float:
        return 0.5 * (1.0 + math.erf(value / math.sqrt(2.0)))

    @classmethod
    def pvalue(cls, value: float) -> float:
        return 2.0 * (1.0 - cls.cdf(abs(value)))

    @staticmethod
    def critical_95() -> float:
        return 1.959963984540054


def _require_pandas():
    try:
        import pandas as pd  # type: ignore
    except ImportError as exc:  # pragma: no cover
        raise RuntimeError("pandas is required for latent geodesic commands.") from exc
    return pd


def _require_numpy():
    try:
        import numpy as np  # type: ignore
    except ImportError as exc:  # pragma: no cover
        raise RuntimeError("numpy is required for latent geodesic commands.") from exc
    return np


def _validate_columns(frame, required_columns: list[str], dataset_name: str) -> None:
    missing = [column for column in required_columns if column not in frame.columns]
    if missing:
        actual = ", ".join(frame.columns)
        needed = ", ".join(missing)
        raise ValueError(f"{dataset_name} is missing required columns: {needed}. Actual columns: {actual}")


def _build_design_matrix_without_distance(frame):
    pd = _require_pandas()
    design = frame[BASE_FEATURE_COLUMNS].copy()
    period_dummies = pd.get_dummies(frame["period_start_year"].astype("string"), prefix="period", dtype=float)
    if not period_dummies.empty:
        period_dummies = period_dummies.drop(columns=[sorted(period_dummies.columns)[0]])
    design = pd.concat([design, period_dummies], axis=1)
    design.insert(0, "const", 1.0)
    return design, list(design.columns) + ["ln_latent_distance_km"]


def _extract_country_names(frame, country_codes: list[str]) -> dict[str, str]:
    names: dict[str, str] = {}
    for code in country_codes:
        origin = frame.loc[frame["origin_iso3"].eq(code), "origin_name"]
        if origin.empty:
            origin = frame.loc[frame["destination_iso3"].eq(code), "destination_name"]
        names[code] = str(origin.iloc[0]) if not origin.empty else code
    return names


def _angles_to_vectors(lat_rad, lon_rad):
    np = _require_numpy()
    cos_lat = np.cos(lat_rad)
    return np.column_stack([cos_lat * np.cos(lon_rad), cos_lat * np.sin(lon_rad), np.sin(lat_rad)])


def _vectors_to_angles(vectors):
    np = _require_numpy()
    safe = vectors / np.clip(np.linalg.norm(vectors, axis=1, keepdims=True), 1e-12, None)
    return np.arcsin(np.clip(safe[:, 2], -1.0, 1.0)), np.arctan2(safe[:, 1], safe[:, 0])


def _build_reference_coordinates(panel, settings: Settings):
    pd = _require_pandas()

    origin = panel[["origin_canonical_id", "origin_iso3", "origin_name"]].rename(columns={"origin_canonical_id": "canonical_id", "origin_iso3": "iso3", "origin_name": "country_name"})
    destination = panel[["destination_canonical_id", "destination_iso3", "destination_name"]].rename(columns={"destination_canonical_id": "canonical_id", "destination_iso3": "iso3", "destination_name": "country_name"})
    countries = pd.concat([origin, destination], ignore_index=True).drop_duplicates(subset=["canonical_id"]).sort_values("iso3").reset_index(drop=True)
    countries["country_index"] = range(len(countries))

    metadata_path = settings.raw_dir / "world_bank" / "country_metadata.csv"
    if not metadata_path.exists():
        raise FileNotFoundError(
            "World Bank country metadata is missing. Run `python -m gravity_world.cli download --source world_bank_wdi` first."
        )
    metadata = pd.read_csv(metadata_path)
    _validate_columns(metadata, ["id", "name", "longitude", "latitude"], "World Bank country metadata")
    metadata = metadata.dropna(subset=["longitude", "latitude"]).copy()
    metadata["longitude"] = pd.to_numeric(metadata["longitude"], errors="coerce")
    metadata["latitude"] = pd.to_numeric(metadata["latitude"], errors="coerce")
    metadata = metadata.dropna(subset=["longitude", "latitude"])
    metadata = metadata.rename(
        columns={
            "id": "iso3",
            "name": "world_bank_name",
            "latitude": "reference_lat_deg",
            "longitude": "reference_lon_deg",
        }
    )
    countries = countries.merge(
        metadata[["iso3", "world_bank_name", "reference_lat_deg", "reference_lon_deg"]],
        how="left",
        on="iso3",
    )
    missing = countries.loc[countries["reference_lat_deg"].isna() | countries["reference_lon_deg"].isna(), "iso3"].tolist()
    if missing:
        raise ValueError(
            "World Bank country metadata is missing coordinates for these panel countries: "
            + ", ".join(sorted(missing))
        )

    return countries


def _angles_to_unconstrained(lat_rad, lon_rad):
    np = _require_numpy()
    return np.arctanh(np.clip(lat_rad / HALF_PI, -0.999999, 0.999999)), np.arctanh(np.clip(lon_rad / math.pi, -0.999999, 0.999999))


def _unconstrained_to_angles(lat_param, lon_param):
    np = _require_numpy()
    lat_share = np.tanh(lat_param)
    lon_share = np.tanh(lon_param)
    return HALF_PI * lat_share, math.pi * lon_share, HALF_PI * (1.0 - lat_share**2), math.pi * (1.0 - lon_share**2)

def _build_initial_beta(settings: Settings, term_names: list[str], mean_flow: float):
    np = _require_numpy()
    beta = np.zeros(len(term_names), dtype=float)
    beta[term_names.index("const")] = math.log(mean_flow if mean_flow > 0 else 1e-12)
    coefficients_path = settings.processed_dir / "models" / "cepii_stock_ppml_coefficients.csv"
    if coefficients_path.exists():
        pd = _require_pandas()
        coefficients = pd.read_csv(coefficients_path)
        coef_map = dict(zip(coefficients["term"], coefficients["estimate"]))
        for idx, term in enumerate(term_names):
            source_term = "ln_distance_km" if term == "ln_latent_distance_km" else term
            if source_term in coef_map:
                beta[idx] = float(coef_map[source_term])
    return beta


def _fit_ppml_coefficients(X, y, beta_init, max_iter: int = 60, tol: float = 1e-9):
    np = _require_numpy()
    beta = beta_init.astype(float).copy()
    eta_clip = 30.0

    def objective(coefficients):
        eta = np.clip(X @ coefficients, -eta_clip, eta_clip)
        mu = np.exp(eta)
        return float((y * eta - mu).sum()), mu

    current_obj, mu = objective(beta)
    converged = False
    iterations = 0
    for iteration in range(1, max_iter + 1):
        xtwx = X.T @ (mu[:, None] * X)
        gradient = X.T @ (y - mu)
        step = np.linalg.pinv(xtwx) @ gradient
        step_scale = 1.0
        accepted = False
        while step_scale >= 1.0 / (2**20):
            candidate_beta = beta + step_scale * step
            candidate_obj, candidate_mu = objective(candidate_beta)
            if candidate_obj >= current_obj - 1e-10:
                beta = candidate_beta
                mu = candidate_mu
                current_obj = candidate_obj
                accepted = True
                break
            step_scale /= 2.0
        if not accepted:
            break
        iterations = iteration
        if float(np.max(np.abs(step_scale * step))) < tol:
            converged = True
            break
    return beta, mu, current_obj, converged, iterations


def _compute_latent_distances(lat_rad, lon_rad, origin_idx, destination_idx):
    np = _require_numpy()
    lat_o = lat_rad[origin_idx]
    lat_d = lat_rad[destination_idx]
    lon_o = lon_rad[origin_idx]
    lon_d = lon_rad[destination_idx]
    delta_lon = lon_o - lon_d
    sin_lat_o = np.sin(lat_o)
    sin_lat_d = np.sin(lat_d)
    cos_lat_o = np.cos(lat_o)
    cos_lat_d = np.cos(lat_d)
    cosine = sin_lat_o * sin_lat_d + cos_lat_o * cos_lat_d * np.cos(delta_lon)
    cosine = np.clip(cosine, -1.0, 1.0)
    central = np.arccos(cosine)
    central_floor = MIN_DISTANCE_KM / EARTH_RADIUS_KM
    active = central > central_floor
    central_safe = np.where(active, central, central_floor)
    distance_km = EARTH_RADIUS_KM * central_safe
    log_distance = np.log(distance_km)
    sin_central = np.sqrt(np.clip(1.0 - cosine**2, 1e-12, None))
    inv_term = np.zeros_like(central_safe)
    inv_term[active] = 1.0 / (central_safe[active] * sin_central[active])
    dcos_dlat_o = cos_lat_o * sin_lat_d - sin_lat_o * cos_lat_d * np.cos(delta_lon)
    dcos_dlat_d = sin_lat_o * cos_lat_d - cos_lat_o * sin_lat_d * np.cos(delta_lon)
    dcos_dlon_o = -cos_lat_o * cos_lat_d * np.sin(delta_lon)
    dcos_dlon_d = cos_lat_o * cos_lat_d * np.sin(delta_lon)
    return (
        log_distance,
        distance_km,
        -inv_term * dcos_dlat_o,
        -inv_term * dcos_dlat_d,
        -inv_term * dcos_dlon_o,
        -inv_term * dcos_dlon_d,
    )


def _compute_quadratic_displacement_penalty(lat_rad, lon_rad, reference_lat_rad, reference_lon_rad, free_indices, penalty_weight: float, observation_count: int):
    np = _require_numpy()
    grad_lat = np.zeros_like(lat_rad)
    grad_lon = np.zeros_like(lon_rad)
    if penalty_weight <= 0.0 or len(free_indices) == 0:
        return 0.0, grad_lat, grad_lon

    current_lat = lat_rad[free_indices]
    current_lon = lon_rad[free_indices]
    reference_lat = reference_lat_rad[free_indices]
    reference_lon = reference_lon_rad[free_indices]
    delta_lon = current_lon - reference_lon
    sin_current_lat = np.sin(current_lat)
    cos_current_lat = np.cos(current_lat)
    sin_reference_lat = np.sin(reference_lat)
    cos_reference_lat = np.cos(reference_lat)
    cosine = sin_current_lat * sin_reference_lat + cos_current_lat * cos_reference_lat * np.cos(delta_lon)
    cosine = np.clip(cosine, -1.0, 1.0)

    scale = penalty_weight * float(observation_count) / float(len(free_indices)) * (EARTH_RADIUS_KM / PENALTY_DISTANCE_SCALE_KM) ** 2
    penalty_value = float(scale * np.sum(2.0 - 2.0 * cosine))

    dcos_dlat = cos_current_lat * sin_reference_lat - sin_current_lat * cos_reference_lat * np.cos(delta_lon)
    dcos_dlon = -cos_current_lat * cos_reference_lat * np.sin(delta_lon)
    grad_lat[free_indices] = 2.0 * scale * dcos_dlat
    grad_lon[free_indices] = 2.0 * scale * dcos_dlon
    return penalty_value, grad_lat, grad_lon


def _build_country_comparison(frame, fitted_total_period_column: str):
    pd = _require_pandas()
    total_years = float(frame[["period_start_year", "period_length_years"]].drop_duplicates()["period_length_years"].sum())
    country_names = _extract_country_names(frame, DEFAULT_COMPARISON_COUNTRIES)
    records = []
    for code in DEFAULT_COMPARISON_COUNTRIES:
        inflow_mask = frame["destination_iso3"].eq(code)
        outflow_mask = frame["origin_iso3"].eq(code)
        observed_inflow_total_period = float(frame.loc[inflow_mask, "flow_total_period"].sum())
        fitted_inflow_total_period = float(frame.loc[inflow_mask, fitted_total_period_column].sum())
        observed_outflow_total_period = float(frame.loc[outflow_mask, "flow_total_period"].sum())
        fitted_outflow_total_period = float(frame.loc[outflow_mask, fitted_total_period_column].sum())
        records.append(
            {
                "country_iso3": code,
                "country_name": country_names.get(code, code),
                "sample_years": total_years,
                "observed_inflow_avg_annual": observed_inflow_total_period / total_years,
                "fitted_inflow_avg_annual": fitted_inflow_total_period / total_years,
                "observed_outflow_avg_annual": observed_outflow_total_period / total_years,
                "fitted_outflow_avg_annual": fitted_outflow_total_period / total_years,
                "observed_balance_avg_annual": (observed_inflow_total_period - observed_outflow_total_period) / total_years,
                "fitted_balance_avg_annual": (fitted_inflow_total_period - fitted_outflow_total_period) / total_years,
            }
        )
    return pd.DataFrame(records)


def _build_spain_importance_label_set(frame, coverage_target: float = 0.95) -> set[str]:
    pd = _require_pandas()
    inbound = frame.loc[frame["destination_iso3"].eq("ESP"), ["origin_iso3", "flow_total_period"]].rename(columns={"origin_iso3": "counterpart_iso3"})
    outbound = frame.loc[frame["origin_iso3"].eq("ESP"), ["destination_iso3", "flow_total_period"]].rename(columns={"destination_iso3": "counterpart_iso3"})
    combined = pd.concat([inbound, outbound], ignore_index=True)
    if combined.empty:
        return {"ESP", "FRA", "USA", "DEU", "ITA", "GBR", "MAR"}
    importance = combined.groupby("counterpart_iso3", as_index=False)["flow_total_period"].sum().sort_values("flow_total_period", ascending=False).reset_index(drop=True)
    total_flow = float(importance["flow_total_period"].sum())
    if total_flow <= 0.0:
        return {"ESP", "FRA", "USA", "DEU", "ITA", "GBR", "MAR"}
    importance["cum_share"] = importance["flow_total_period"].cumsum() / total_flow
    selected = importance.loc[importance["cum_share"].le(coverage_target), "counterpart_iso3"].tolist()
    if len(selected) < len(importance):
        cutoff_index = min(len(selected), len(importance) - 1)
        selected.append(str(importance.iloc[cutoff_index]["counterpart_iso3"]))
    return set(selected) | {"ESP", "FRA", "USA", "DEU", "ITA", "GBR", "MAR"}


def _write_markdown_summary(path: Path, fit_stats: dict, coefficients, country_comparison, spain_by_period) -> None:
    title = fit_stats.get("model_title", "Latent Geodesic PPML")
    lines = [
        f"# {title}",
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
        f"| Observations used in estimation | {fit_stats['n_obs_estimation_sample']:,} |",
        f"| Zero-flow observations included | {fit_stats['n_obs_zero_flow_included']:,} |",
        f"| Parameters | {fit_stats['n_parameters']} |",
        f"| Coordinate iterations | {fit_stats['coordinate_iterations']} |",
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
        "## Coefficients",
        "",
        "| Term | Estimate | Robust Std. Error | z-stat | p-value | 95% CI Low | 95% CI High |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ])
    for row in coefficients.itertuples(index=False):
        lines.append(f"| {row.term} | {row.estimate:.6f} | {row.std_error:.6f} | {row.z_stat:.6f} | {row.p_value_normal_approx:.6g} | {row.ci95_lower:.6f} | {row.ci95_upper:.6f} |")
    lines.extend([
        "",
        "## Average Annual Flow Comparison Over Full Sample",
        "",
        "Fitted values below are latent-distance PPML fitted means.",
        "",
        "| Country | Observed Inflow | Fitted Inflow Mean | Observed Outflow | Fitted Outflow Mean | Observed Balance | Fitted Balance Mean |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ])
    for row in country_comparison.itertuples(index=False):
        lines.append(f"| {row.country_name} ({row.country_iso3}) | {row.observed_inflow_avg_annual:,.3f} | {row.fitted_inflow_avg_annual:,.3f} | {row.observed_outflow_avg_annual:,.3f} | {row.fitted_outflow_avg_annual:,.3f} | {row.observed_balance_avg_annual:,.3f} | {row.fitted_balance_avg_annual:,.3f} |")
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
        lines.append(f"| {period_label} | {row.observed_inflow_spain:,.3f} | {row.fitted_inflow_spain:,.3f} | {row.observed_outflow_spain:,.3f} | {row.fitted_outflow_spain:,.3f} | {observed_balance:,.3f} | {fitted_balance:,.3f} |")
    path.write_text("\n".join(lines), encoding="utf-8")


def _svg_escape(text: str) -> str:
    return str(text).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;").replace('"', "&quot;")


def _write_map_svg(path: Path, coordinates, anchor_codes: tuple[str, ...], spain_label_codes: set[str] | None = None) -> None:
    width = 1900
    height = 980
    left = 60
    top = 50
    plot_width = width - 2 * left
    plot_height = height - 2 * top

    def x_scale(lon: float) -> float:
        return left + (lon + 180.0) / 360.0 * plot_width

    def y_scale(lat: float) -> float:
        return top + (90.0 - lat) / 180.0 * plot_height

    ranked = coordinates.sort_values("displacement_km", ascending=False).reset_index(drop=True)
    label_set = set(ranked.head(20)["iso3"]) | set(spain_label_codes or set()) | set(anchor_codes) | {"ESP", "USA", "DEU", "ITA", "GBR", "MAR"}
    parts = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">', '<rect width="100%" height="100%" fill="#f6f4ed" />', '<style>text { font-family: Segoe UI, Arial, sans-serif; fill: #1f2933; } .title { font-size: 28px; font-weight: 700; } .subtitle { font-size: 14px; fill: #52606d; } .grid { stroke: #d9d5c9; stroke-width: 1; } .ref { fill: #b8bfc6; opacity: 0.8; } .latent { fill: #d94841; opacity: 0.92; } .anchor { fill: #0f4c81; } .label { font-size: 9px; } .legend { font-size: 12px; }</style>']
    parts.append(f'<text x="{left}" y="30" class="title">Latent Geodesic Migration Map</text>')
    anchor_text = ', '.join(anchor_codes)
    parts.append(f'<text x="{left}" y="48" class="subtitle">Gray = World Bank capital coordinates. Red = latent PPML coordinates. Fixed anchors: {anchor_text}.</text>')
    for lon in range(-150, 181, 30):
        x = x_scale(float(lon))
        parts.append(f'<line x1="{x:.2f}" y1="{top}" x2="{x:.2f}" y2="{top + plot_height}" class="grid" />')
        parts.append(f'<text x="{x:.2f}" y="{top + plot_height + 18}" text-anchor="middle" class="legend">{lon}</text>')
    for lat in range(-60, 91, 30):
        y = y_scale(float(lat))
        parts.append(f'<line x1="{left}" y1="{y:.2f}" x2="{left + plot_width}" y2="{y:.2f}" class="grid" />')
        parts.append(f'<text x="{left - 8}" y="{y + 4:.2f}" text-anchor="end" class="legend">{lat}</text>')
    for row in coordinates.itertuples(index=False):
        x0 = x_scale(float(row.reference_lon_deg)); y0 = y_scale(float(row.reference_lat_deg)); x1 = x_scale(float(row.latent_lon_deg)); y1 = y_scale(float(row.latent_lat_deg))
        stroke = "#c1c7cd" if not bool(row.is_anchor) else "#7b8794"
        width_px = 1.0 if not bool(row.is_anchor) else 1.8
        parts.append(f'<line x1="{x0:.2f}" y1="{y0:.2f}" x2="{x1:.2f}" y2="{y1:.2f}" stroke="{stroke}" stroke-width="{width_px}" opacity="0.7" />')
        ref_class = "anchor" if bool(row.is_anchor) else "ref"
        latent_class = "anchor" if bool(row.is_anchor) else "latent"
        parts.append(f'<circle cx="{x0:.2f}" cy="{y0:.2f}" r="2.2" class="{ref_class}" />')
        parts.append(f'<circle cx="{x1:.2f}" cy="{y1:.2f}" r="2.8" class="{latent_class}" />')
        if row.iso3 in label_set:
            parts.append(f'<text x="{x1 + 4:.2f}" y="{y1 - 4:.2f}" class="label">{_svg_escape(row.iso3)}</text>')
    legend_x = left + 20
    legend_y = top + 18
    parts.append(f'<circle cx="{legend_x}" cy="{legend_y}" r="3" class="ref" />')
    parts.append(f'<text x="{legend_x + 12}" y="{legend_y + 4}" class="legend">Reference position</text>')
    parts.append(f'<circle cx="{legend_x + 160}" cy="{legend_y}" r="3" class="latent" />')
    parts.append(f'<text x="{legend_x + 172}" y="{legend_y + 4}" class="legend">Latent PPML position</text>')
    parts.append('</svg>')
    path.write_text("\n".join(parts), encoding="utf-8")

def estimate_latent_geodesic_ppml(settings: Settings, output_dir: Path | None = None, penalty_weight: float = 0.0, anchors: tuple[str, ...] | list[str] | None = None) -> list[Path]:
    pd = _require_pandas()
    np = _require_numpy()
    if penalty_weight < 0.0:
        raise ValueError("penalty_weight must be nonnegative.")

    anchor_codes = tuple(dict.fromkeys(code.upper() for code in (anchors or ("ESP", "FRA"))))
    if len(anchor_codes) < 2:
        raise ValueError("At least two anchor countries are required.")

    penalized = penalty_weight > 0.0
    model_prefix = "latent_geodesic_penalized_ppml" if penalized else "latent_geodesic_ppml"
    curated_base = "latent_geodesic_penalized" if penalized else "latent_geodesic"
    summary_title = "Penalized Latent Geodesic PPML" if penalized else "Latent Geodesic PPML"

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
        raise ValueError("No observations available for latent geodesic PPML after filtering.")

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

    design_base, term_names = _build_design_matrix_without_distance(estimation)
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

    lat_param, lon_param = _angles_to_unconstrained(reference_lat_rad[free_indices], reference_lon_rad[free_indices])
    beta = _build_initial_beta(settings, term_names, float(y.mean()))
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

    reference_log_distance, _, _, _, _, _ = _compute_latent_distances(reference_lat_rad, reference_lon_rad, origin_idx, destination_idx)
    reference_X = np.column_stack([X_base, reference_log_distance])
    _, _, reference_obj_no_const, _, _ = _fit_ppml_coefficients(reference_X, y, beta.copy())

    for iteration in range(1, 61):
        lat_rad, lon_rad, dlat_param, dlon_param = compose_angles(lat_param, lon_param)
        log_distance, distance_km, dlat_o, dlat_d, dlon_o, dlon_d = _compute_latent_distances(lat_rad, lon_rad, origin_idx, destination_idx)
        X = np.column_stack([X_base, log_distance])
        beta, mu, objective_no_const, beta_converged, beta_iterations = _fit_ppml_coefficients(X, y, beta)
        penalty_value, penalty_grad_lat, penalty_grad_lon = _compute_quadratic_displacement_penalty(
            lat_rad, lon_rad, reference_lat_rad, reference_lon_rad, free_indices, penalty_weight, len(y)
        )
        penalized_objective = objective_no_const - penalty_value
        objective_history.append(penalized_objective)
        if best_state is None or penalized_objective > float(best_state["objective"]):
            best_state = {
                "objective": float(penalized_objective),
                "objective_no_const": float(objective_no_const),
                "penalty_value": float(penalty_value),
                "beta": beta.copy(),
                "mu": mu.copy(),
                "log_distance": log_distance.copy(),
                "distance_km": distance_km.copy(),
                "lat_rad": lat_rad.copy(),
                "lon_rad": lon_rad.copy(),
                "iteration": iteration,
                "beta_converged": beta_converged,
                "beta_iterations": beta_iterations,
            }
        grad_common = (y - mu) * float(beta[-1])
        grad_lat = np.bincount(origin_idx, weights=grad_common * dlat_o, minlength=country_count) + np.bincount(destination_idx, weights=grad_common * dlat_d, minlength=country_count)
        grad_lon = np.bincount(origin_idx, weights=grad_common * dlon_o, minlength=country_count) + np.bincount(destination_idx, weights=grad_common * dlon_d, minlength=country_count)
        grad_lat_total = grad_lat + penalty_grad_lat
        grad_lon_total = grad_lon + penalty_grad_lon
        grad_lat_param = np.clip(grad_lat_total[free_indices] * dlat_param, -5000.0, 5000.0)
        grad_lon_param = np.clip(grad_lon_total[free_indices] * dlon_param, -5000.0, 5000.0)
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
            break

    if best_state is None:
        raise RuntimeError("Latent geodesic PPML failed to produce any evaluated state.")

    final_beta = np.asarray(best_state["beta"], dtype=float)
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
    z_stats = final_beta / std_errors
    critical = _NormalApprox.critical_95()
    p_values = np.array([_NormalApprox.pvalue(float(value)) for value in z_stats])
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
    estimation["fitted_flow"] = final_mu
    estimation["fitted_flow_total_period"] = final_mu * period_lengths
    estimation["raw_residual_flow"] = y - final_mu

    coefficients = pd.DataFrame({
        "term": term_names,
        "estimate": final_beta,
        "std_error": std_errors,
        "z_stat": z_stats,
        "p_value_normal_approx": p_values,
        "ci95_lower": ci_lower,
        "ci95_upper": ci_upper,
    })
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
        "prediction_target": "E[flow|X] = exp(Xb)",
        "n_countries": int(country_count),
        "n_free_countries": int(len(free_indices)),
        "anchors": ", ".join(anchor_codes),
        "n_obs_total_panel": int(len(panel)),
        "n_obs_estimation_sample": int(len(estimation)),
        "n_obs_zero_flow_included": int((estimation["flow"] == 0).sum()),
        "n_parameters": int(len(final_beta) + 2 * len(free_indices)),
        "coordinate_iterations": final_iteration,
        "beta_converged_last_step": bool(best_state["beta_converged"]),
        "beta_iterations_last_step": int(best_state["beta_iterations"]),
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
    fitted_sample_path = output_dir / f"{model_prefix}_fitted_sample.csv"
    spain_by_period_path = output_dir / f"{model_prefix}_spain_by_period.csv"
    country_comparison_path = output_dir / f"{model_prefix}_country_comparison.csv"
    coordinates_path = output_dir / f"{model_prefix}_country_coordinates.csv"
    map_path = output_dir / f"{model_prefix}_map.svg"

    summary_lines = [
        f"{summary_title} Summary",
        "",
        f"Observations used in estimation: {len(estimation):,}",
        f"Countries: {country_count}",
        f"Free countries: {len(free_indices)}",
        f"Anchors: {fit_stats['anchors']}",
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
        summary_lines.extend([
            f"Quadratic penalty weight: {penalty_weight:.6f}",
            f"Quadratic penalty value: {final_penalty_value:.3f}",
            f"Penalized objective: {penalized_objective:.3f}",
            f"Penalized objective improvement: {penalized_objective - reference_log_likelihood:.3f}",
        ])
    summary_lines.extend([
        "",
        "Coefficients:",
        "term,estimate,std_error,z_stat,p_value_normal_approx,ci95_lower,ci95_upper",
        *[f"{row.term},{row.estimate:.6f},{row.std_error:.6f},{row.z_stat:.6f},{row.p_value_normal_approx:.6g},{row.ci95_lower:.6f},{row.ci95_upper:.6f}" for row in coefficients.itertuples(index=False)],
    ])
    summary_txt_path.write_text("\n".join(summary_lines), encoding="utf-8")
    fit_stats_path.write_text(json.dumps(fit_stats, indent=2), encoding="utf-8")
    coefficients.to_csv(coefficients_path, index=False)
    estimation.to_csv(fitted_sample_path, index=False)
    spain_by_period.to_csv(spain_by_period_path, index=False)
    country_comparison.to_csv(country_comparison_path, index=False)
    coordinates.to_csv(coordinates_path, index=False)
    _write_markdown_summary(summary_md_path, fit_stats, coefficients, country_comparison, spain_by_period)
    _write_map_svg(map_path, coordinates, anchor_codes=anchor_codes, spain_label_codes=spain_label_codes)

    curated_dir = settings.root_dir / "output"
    curated_dir.mkdir(parents=True, exist_ok=True)
    curated_summary_path = curated_dir / f"{curated_base}_summary.md"
    curated_map_path = curated_dir / f"{curated_base}_map.svg"
    shutil.copyfile(summary_md_path, curated_summary_path)
    shutil.copyfile(map_path, curated_map_path)
    return [
        summary_txt_path, summary_md_path, fit_stats_path, coefficients_path, fitted_sample_path,
        spain_by_period_path, country_comparison_path, coordinates_path, map_path, curated_summary_path, curated_map_path,
    ]
