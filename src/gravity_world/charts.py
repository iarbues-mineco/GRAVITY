from __future__ import annotations

import math
from pathlib import Path

from .settings import Settings


SVG_WIDTH = 1600
LEFT_MARGIN = 260
RIGHT_MARGIN = 40
TOP_MARGIN = 80
BOTTOM_MARGIN = 40
ROW_HEIGHT = 20
BAR_HEIGHT = 7
BAR_GAP = 2
TITLE_HEIGHT = 28
SUBTITLE_HEIGHT = 18


def _require_pandas():
    try:
        import pandas as pd  # type: ignore
    except ImportError as exc:  # pragma: no cover
        raise RuntimeError("pandas is required for chart commands.") from exc
    return pd


def _svg_escape(text: str) -> str:
    return (
        str(text)
        .replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace('"', "&quot;")
    )


def _load_fitted_sample(settings: Settings, model_prefix: str):
    pd = _require_pandas()
    path = settings.processed_dir / "models" / f"{model_prefix}_fitted_sample.csv"
    if not path.exists():
        raise FileNotFoundError(
            f"Fitted sample for `{model_prefix}` is missing. Run the corresponding model estimation command first."
        )
    return pd.read_csv(path)


def build_all_country_inflow_comparison(
    settings: Settings,
    model_prefix: str = "cepii_gravity_ols",
    scale: str = "linear",
) -> list[Path]:
    pd = _require_pandas()

    if scale not in {"linear", "log"}:
        raise ValueError("scale must be either `linear` or `log`.")

    fitted = _load_fitted_sample(settings, model_prefix)
    required = ["destination_iso3", "destination_name", "flow_total_period", "fitted_flow_total_period", "period_length_years"]
    missing = [column for column in required if column not in fitted.columns]
    if missing:
        actual = ", ".join(fitted.columns)
        needed = ", ".join(missing)
        raise ValueError(f"Fitted sample is missing required columns: {needed}. Actual columns: {actual}")

    total_years = float(fitted[["period_start_year", "period_length_years"]].drop_duplicates()["period_length_years"].sum())
    inflows = (
        fitted.groupby(["destination_iso3", "destination_name"], as_index=False)[["flow_total_period", "fitted_flow_total_period"]]
        .sum()
        .rename(columns={
            "destination_iso3": "country_iso3",
            "destination_name": "country_name",
            "flow_total_period": "observed_inflow_total_period",
            "fitted_flow_total_period": "fitted_inflow_total_period",
        })
    )
    inflows["sample_years"] = total_years
    inflows["observed_inflow_avg_annual"] = inflows["observed_inflow_total_period"] / total_years
    inflows["fitted_inflow_avg_annual"] = inflows["fitted_inflow_total_period"] / total_years
    inflows["fitted_to_observed_ratio"] = inflows["fitted_inflow_avg_annual"] / inflows["observed_inflow_avg_annual"].replace({0: pd.NA})
    inflows = inflows.sort_values(["observed_inflow_avg_annual", "country_name"], ascending=[False, True]).reset_index(drop=True)

    output_dir = settings.processed_dir / "models"
    output_dir.mkdir(parents=True, exist_ok=True)
    csv_path = output_dir / f"{model_prefix}_all_country_inflow_comparison.csv"
    suffix = "" if scale == "linear" else "_log"
    svg_path = output_dir / f"{model_prefix}_all_country_inflow_comparison{suffix}.svg"

    inflows.to_csv(csv_path, index=False)
    _write_svg_chart(svg_path, inflows, model_prefix, total_years, scale=scale)
    return [csv_path, svg_path]


def _write_svg_chart(path: Path, inflows, model_prefix: str, total_years: float, scale: str = "linear") -> None:
    n = len(inflows)
    chart_width = SVG_WIDTH - LEFT_MARGIN - RIGHT_MARGIN
    chart_height = n * ROW_HEIGHT
    height = TOP_MARGIN + chart_height + BOTTOM_MARGIN
    max_value = float(max(inflows["observed_inflow_avg_annual"].max(), inflows["fitted_inflow_avg_annual"].max(), 1.0))
    x_ticks = 5
    use_log_scale = scale == "log"

    def scale_value(value: float) -> float:
        if max_value <= 0:
            return 0.0
        if use_log_scale:
            return math.log10(value + 1.0) / math.log10(max_value + 1.0)
        return value / max_value

    parts: list[str] = []
    parts.append(f'<svg xmlns="http://www.w3.org/2000/svg" width="{SVG_WIDTH}" height="{height}" viewBox="0 0 {SVG_WIDTH} {height}">')
    parts.append('<rect width="100%" height="100%" fill="#fffdf7" />')
    parts.append('<style>')
    parts.append('text { font-family: Segoe UI, Arial, sans-serif; fill: #1f2933; }')
    parts.append('.title { font-size: 24px; font-weight: 700; }')
    parts.append('.subtitle { font-size: 13px; fill: #52606d; }')
    parts.append('.label { font-size: 11px; }')
    parts.append('.tick { font-size: 10px; fill: #52606d; }')
    parts.append('.legend { font-size: 11px; }')
    parts.append('</style>')

    parts.append(f'<text x="{LEFT_MARGIN}" y="32" class="title">Observed vs Fitted Inflows by Country</text>')
    parts.append(
        f'<text x="{LEFT_MARGIN}" y="52" class="subtitle">Model: {_svg_escape(model_prefix)} | Metric: average annual inflow over {int(total_years)} years | Sorted by observed inflow | Scale: {_svg_escape(scale)}</text>'
    )
    parts.append(f'<rect x="{LEFT_MARGIN}" y="60" width="12" height="12" fill="#1f77b4" rx="2" ry="2" />')
    parts.append(f'<text x="{LEFT_MARGIN + 18}" y="70" class="legend">Observed</text>')
    parts.append(f'<rect x="{LEFT_MARGIN + 90}" y="60" width="12" height="12" fill="#e76f51" rx="2" ry="2" />')
    parts.append(f'<text x="{LEFT_MARGIN + 108}" y="70" class="legend">Fitted</text>')

    for tick in range(x_ticks + 1):
        if use_log_scale:
            share = tick / x_ticks
            value = (10 ** (math.log10(max_value + 1.0) * share)) - 1.0
        else:
            value = max_value * tick / x_ticks
        x = LEFT_MARGIN + chart_width * tick / x_ticks
        parts.append(f'<line x1="{x:.2f}" y1="{TOP_MARGIN}" x2="{x:.2f}" y2="{TOP_MARGIN + chart_height}" stroke="#e5e7eb" stroke-width="1" />')
        parts.append(f'<text x="{x:.2f}" y="{TOP_MARGIN - 8}" text-anchor="middle" class="tick">{value:,.0f}</text>')

    for idx, row in inflows.iterrows():
        y = TOP_MARGIN + idx * ROW_HEIGHT
        observed = float(row["observed_inflow_avg_annual"])
        fitted = float(row["fitted_inflow_avg_annual"])
        observed_width = chart_width * scale_value(observed)
        fitted_width = chart_width * scale_value(fitted)
        label = f"{row['country_name']} ({row['country_iso3']})"
        parts.append(f'<text x="{LEFT_MARGIN - 8}" y="{y + 13}" text-anchor="end" class="label">{_svg_escape(label)}</text>')
        parts.append(f'<rect x="{LEFT_MARGIN}" y="{y + 3}" width="{observed_width:.2f}" height="{BAR_HEIGHT}" fill="#1f77b4" />')
        parts.append(f'<rect x="{LEFT_MARGIN}" y="{y + 3 + BAR_HEIGHT + BAR_GAP}" width="{fitted_width:.2f}" height="{BAR_HEIGHT}" fill="#e76f51" />')

    parts.append('</svg>')
    path.write_text("\n".join(parts), encoding="utf-8")
