from __future__ import annotations

import json
import math
import shutil
import struct
from pathlib import Path

from .settings import Settings

DEFAULT_LATENT_MODEL_PREFIXES = ["latent_geodesic_ppml", "latent_geodesic_penalized_ppml"]
DEFAULT_FIXED_LABELS = {"ESP", "USA", "DEU", "ITA", "GBR", "MAR", "FRA", "CHN"}


def _require_pandas():
    try:
        import pandas as pd  # type: ignore
    except ImportError as exc:  # pragma: no cover
        raise RuntimeError("pandas is required for cartogram commands.") from exc
    return pd


def _strip_dbf_text(value: bytes) -> str:
    return value.decode("latin1").replace("\x00", "").strip()


def _read_dbf_records(path: Path) -> list[dict[str, str]]:
    with path.open("rb") as handle:
        header = handle.read(32)
        record_count = struct.unpack("<I", header[4:8])[0]
        header_length = struct.unpack("<H", header[8:10])[0]
        record_length = struct.unpack("<H", header[10:12])[0]
        fields: list[tuple[str, str, int]] = []
        while True:
            descriptor = handle.read(32)
            if not descriptor or descriptor[0] == 0x0D:
                break
            name = descriptor[:11].split(b"\x00", 1)[0].decode("ascii")
            field_type = chr(descriptor[11])
            field_length = descriptor[16]
            fields.append((name, field_type, field_length))
        handle.seek(header_length)
        rows: list[dict[str, str]] = []
        for _ in range(record_count):
            raw_record = handle.read(record_length)
            if not raw_record or raw_record[0] == 0x2A:
                continue
            offset = 1
            row: dict[str, str] = {}
            for name, _field_type, field_length in fields:
                row[name] = _strip_dbf_text(raw_record[offset : offset + field_length])
                offset += field_length
            rows.append(row)
    return rows


def _read_polygon_shapes(path: Path) -> list[list[list[tuple[float, float]]]]:
    shapes: list[list[list[tuple[float, float]]]] = []
    with path.open("rb") as handle:
        handle.read(100)
        while True:
            record_header = handle.read(8)
            if not record_header:
                break
            if len(record_header) != 8:
                raise ValueError("Unexpected truncated shapefile record header.")
            _record_number, content_length_words = struct.unpack(">2i", record_header)
            content = handle.read(content_length_words * 2)
            if len(content) != content_length_words * 2:
                raise ValueError("Unexpected truncated shapefile record body.")
            shape_type = struct.unpack("<i", content[:4])[0]
            if shape_type == 0:
                shapes.append([])
                continue
            if shape_type != 5:
                raise ValueError(f"Unsupported shapefile shape type: {shape_type}")
            num_parts = struct.unpack("<i", content[36:40])[0]
            num_points = struct.unpack("<i", content[40:44])[0]
            part_indices = list(struct.unpack(f"<{num_parts}i", content[44 : 44 + 4 * num_parts]))
            points_offset = 44 + 4 * num_parts
            points: list[tuple[float, float]] = []
            for point_idx in range(num_points):
                point_offset = points_offset + 16 * point_idx
                lon, lat = struct.unpack("<2d", content[point_offset : point_offset + 16])
                points.append((float(lon), float(lat)))
            rings: list[list[tuple[float, float]]] = []
            part_ends = part_indices[1:] + [num_points]
            for start, end in zip(part_indices, part_ends):
                ring = points[start:end]
                if ring:
                    rings.append(ring)
            shapes.append(rings)
    return shapes


def _load_country_geometries(settings: Settings, needed_codes: set[str]) -> dict[str, list[list[tuple[float, float]]]]:
    base = settings.raw_dir / "maps" / "ne_50m_admin_0_countries"
    shp_path = base / "ne_50m_admin_0_countries.shp"
    dbf_path = base / "ne_50m_admin_0_countries.dbf"
    if not shp_path.exists() or not dbf_path.exists():
        raise FileNotFoundError(
            "Natural Earth country polygons are missing. Expected files under data/raw/maps/ne_50m_admin_0_countries/."
        )
    records = _read_dbf_records(dbf_path)
    shapes = _read_polygon_shapes(shp_path)
    if len(records) != len(shapes):
        raise ValueError("Shapefile geometry count does not match DBF record count.")

    geometries: dict[str, list[list[tuple[float, float]]]] = {}
    for record, rings in zip(records, shapes):
        code = record.get("ADM0_A3", "")
        if code in needed_codes:
            geometries.setdefault(code, []).extend(rings)
    missing = sorted(needed_codes - set(geometries))
    if missing:
        raise ValueError("Natural Earth polygons are missing these model countries: " + ", ".join(missing))
    return geometries


def _ring_area_and_centroid(ring: list[tuple[float, float]]) -> tuple[float, float, float]:
    if len(ring) < 3:
        if not ring:
            return 0.0, 0.0, 0.0
        mean_lon = sum(point[0] for point in ring) / len(ring)
        mean_lat = sum(point[1] for point in ring) / len(ring)
        return 0.0, mean_lon, mean_lat
    closed = ring if ring[0] == ring[-1] else ring + [ring[0]]
    cross_sum = 0.0
    centroid_x = 0.0
    centroid_y = 0.0
    for (x0, y0), (x1, y1) in zip(closed[:-1], closed[1:]):
        cross = x0 * y1 - x1 * y0
        cross_sum += cross
        centroid_x += (x0 + x1) * cross
        centroid_y += (y0 + y1) * cross
    if abs(cross_sum) < 1e-12:
        mean_lon = sum(point[0] for point in ring) / len(ring)
        mean_lat = sum(point[1] for point in ring) / len(ring)
        return 0.0, mean_lon, mean_lat
    area = 0.5 * cross_sum
    return area, centroid_x / (3.0 * cross_sum), centroid_y / (3.0 * cross_sum)


def _geometry_centroid(rings: list[list[tuple[float, float]]]) -> tuple[float, float]:
    area_total = 0.0
    centroid_lon_total = 0.0
    centroid_lat_total = 0.0
    all_points: list[tuple[float, float]] = []
    for ring in rings:
        all_points.extend(ring)
        area, centroid_lon, centroid_lat = _ring_area_and_centroid(ring)
        area_total += area
        centroid_lon_total += centroid_lon * area
        centroid_lat_total += centroid_lat * area
    if abs(area_total) < 1e-12:
        if not all_points:
            return 0.0, 0.0
        return (
            sum(point[0] for point in all_points) / len(all_points),
            sum(point[1] for point in all_points) / len(all_points),
        )
    return centroid_lon_total / area_total, centroid_lat_total / area_total


def _build_spain_label_set(settings: Settings, coverage_target: float = 0.95) -> set[str]:
    pd = _require_pandas()
    panel_path = settings.processed_dir / "panels" / "bilateral_panel_cepii_stock.csv"
    panel = pd.read_csv(panel_path, usecols=["origin_iso3", "destination_iso3", "flow_total_period"])
    inbound = panel.loc[panel["destination_iso3"].eq("ESP"), ["origin_iso3", "flow_total_period"]].rename(columns={"origin_iso3": "counterpart_iso3"})
    outbound = panel.loc[panel["origin_iso3"].eq("ESP"), ["destination_iso3", "flow_total_period"]].rename(columns={"destination_iso3": "counterpart_iso3"})
    combined = pd.concat([inbound, outbound], ignore_index=True)
    if combined.empty:
        return set(DEFAULT_FIXED_LABELS)
    importance = combined.groupby("counterpart_iso3", as_index=False)["flow_total_period"].sum().sort_values("flow_total_period", ascending=False).reset_index(drop=True)
    total_flow = float(importance["flow_total_period"].sum())
    if total_flow <= 0.0:
        return set(DEFAULT_FIXED_LABELS)
    importance["cum_share"] = importance["flow_total_period"].cumsum() / total_flow
    selected = importance.loc[importance["cum_share"].le(coverage_target), "counterpart_iso3"].tolist()
    if len(selected) < len(importance):
        cutoff_index = min(len(selected), len(importance) - 1)
        selected.append(str(importance.iloc[cutoff_index]["counterpart_iso3"]))
    return set(selected) | set(DEFAULT_FIXED_LABELS)


def _compute_bounds(records: list[dict[str, object]]) -> tuple[float, float, float, float]:
    lon_values: list[float] = []
    lat_values: list[float] = []
    for record in records:
        delta_lon = float(record["delta_lon"])
        delta_lat = float(record["delta_lat"])
        centroid_lon = float(record["centroid_lon"])
        centroid_lat = float(record["centroid_lat"])
        lon_values.extend([centroid_lon, centroid_lon + delta_lon])
        lat_values.extend([centroid_lat, centroid_lat + delta_lat])
        for ring in record["rings"]:  # type: ignore[index]
            for lon, lat in ring:
                lon_values.append(float(lon))
                lat_values.append(float(lat))
                lon_values.append(float(lon + delta_lon))
                lat_values.append(float(lat + delta_lat))
    lon_min = min(lon_values)
    lon_max = max(lon_values)
    lat_min = min(lat_values)
    lat_max = max(lat_values)
    lon_span = max(lon_max - lon_min, 1.0)
    lat_span = max(lat_max - lat_min, 1.0)
    lon_pad = max(10.0, 0.04 * lon_span)
    lat_pad = max(6.0, 0.05 * lat_span)
    return lon_min - lon_pad, lon_max + lon_pad, lat_min - lat_pad, lat_max + lat_pad


def _project_factory(bounds: tuple[float, float, float, float], width: int, height: int, left: int, top: int):
    lon_min, lon_max, lat_min, lat_max = bounds
    plot_width = width - 2 * left
    plot_height = height - 2 * top
    def project(lon: float, lat: float) -> tuple[float, float]:
        x = left + (lon - lon_min) / (lon_max - lon_min) * plot_width
        y = top + (lat_max - lat) / (lat_max - lat_min) * plot_height
        return x, y
    return project


def _rings_to_path(rings: list[list[tuple[float, float]]], project, delta_lon: float = 0.0, delta_lat: float = 0.0) -> str:
    commands: list[str] = []
    for ring in rings:
        if not ring:
            continue
        first_x, first_y = project(ring[0][0] + delta_lon, ring[0][1] + delta_lat)
        commands.append(f"M {first_x:.2f} {first_y:.2f}")
        for lon, lat in ring[1:]:
            x, y = project(lon + delta_lon, lat + delta_lat)
            commands.append(f"L {x:.2f} {y:.2f}")
        commands.append("Z")
    return " ".join(commands)


def _svg_escape(text: str) -> str:
    return str(text).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;").replace('"', "&quot;")


def _load_map_records(settings: Settings, model_prefix: str) -> tuple[list[dict[str, object]], tuple[str, ...], str]:
    pd = _require_pandas()
    coordinates_path = settings.processed_dir / "models" / f"{model_prefix}_country_coordinates.csv"
    if not coordinates_path.exists():
        raise FileNotFoundError(f"Country-coordinate file is missing for model prefix `{model_prefix}`.")
    coordinates = pd.read_csv(coordinates_path)
    needed_columns = {
        "iso3", "country_name", "reference_lat_deg", "reference_lon_deg", "latent_lat_deg", "latent_lon_deg", "is_anchor"
    }
    missing_columns = sorted(needed_columns - set(coordinates.columns))
    if missing_columns:
        raise ValueError("Coordinate file is missing required columns: " + ", ".join(missing_columns))

    important_labels = _build_spain_label_set(settings)
    anchor_codes = tuple(coordinates.loc[coordinates["is_anchor"].astype(bool), "iso3"].astype(str).tolist())
    geometries = _load_country_geometries(settings, set(coordinates["iso3"].astype(str)))

    records: list[dict[str, object]] = []
    for row in coordinates.itertuples(index=False):
        iso3 = str(row.iso3)
        rings = geometries[iso3]
        centroid_lon, centroid_lat = _geometry_centroid(rings)
        delta_lon = float(row.latent_lon_deg) - float(row.reference_lon_deg)
        delta_lat = float(row.latent_lat_deg) - float(row.reference_lat_deg)
        records.append(
            {
                "iso3": iso3,
                "country_name": str(row.country_name),
                "rings": rings,
                "centroid_lon": centroid_lon,
                "centroid_lat": centroid_lat,
                "delta_lon": delta_lon,
                "delta_lat": delta_lat,
                "target_lon": centroid_lon + delta_lon,
                "target_lat": centroid_lat + delta_lat,
                "is_anchor": bool(row.is_anchor),
                "is_spain_relevant": iso3 in important_labels,
            }
        )
    return records, anchor_codes, model_prefix.replace("_ppml", "")


def _fill_color(record: dict[str, object], shifted: bool) -> str:
    if bool(record["is_anchor"]):
        return "#d9e8f5" if not shifted else "#a9c9e4"
    if bool(record["is_spain_relevant"]):
        return "#e4f0db" if not shifted else "#f3cdbd"
    return "#edf0e8" if not shifted else "#f7ddd2"


def _stroke_color(record: dict[str, object], shifted: bool) -> str:
    if bool(record["is_anchor"]):
        return "#3d6f98"
    return "#aeb5aa" if not shifted else "#bf7a60"


def _arrow_color(record: dict[str, object]) -> str:
    if bool(record["is_anchor"]):
        return "#3d6f98"
    if bool(record["is_spain_relevant"]):
        return "#8c4545"
    return "#8c9298"


def _draw_legend(parts: list[str], left: int, top: int, shifted: bool) -> None:
    legend_x = left + 8
    legend_y = top + 12
    parts.append(f'<rect x="{legend_x - 4}" y="{legend_y - 16}" width="350" height="62" rx="8" fill="#ffffff" opacity="0.78" />')
    parts.append(f'<rect x="{legend_x}" y="{legend_y}" width="18" height="10" fill="#edf0e8" stroke="#aeb5aa" />')
    parts.append(f'<text x="{legend_x + 26}" y="{legend_y + 9}" class="legend">Real country surface</text>')
    if shifted:
        parts.append(f'<rect x="{legend_x + 150}" y="{legend_y}" width="18" height="10" fill="#f7ddd2" stroke="#bf7a60" />')
        parts.append(f'<text x="{legend_x + 176}" y="{legend_y + 9}" class="legend">Shifted country surface</text>')
    parts.append(f'<line x1="{legend_x}" y1="{legend_y + 26}" x2="{legend_x + 26}" y2="{legend_y + 26}" stroke="#8c4545" stroke-width="1.8" marker-end="url(#arrowhead)" />')
    parts.append(f'<text x="{legend_x + 34}" y="{legend_y + 30}" class="legend">Centroid displacement arrow</text>')


def _write_static_arrows_map(path: Path, records: list[dict[str, object]], anchor_codes: tuple[str, ...], model_label: str) -> None:
    width = 1960
    height = 1080
    left = 40
    top = 60
    bounds = _compute_bounds(records)
    project = _project_factory(bounds, width, height, left, top)
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#f7f3ea" />',
        '<defs><marker id="arrowhead" markerWidth="8" markerHeight="8" refX="7" refY="3.5" orient="auto"><polygon points="0 0, 8 3.5, 0 7" fill="#8c4545" /></marker></defs>',
        '<style>text { font-family: Segoe UI, Arial, sans-serif; fill: #243039; } .title { font-size: 28px; font-weight: 700; } .subtitle { font-size: 14px; fill: #57636c; } .label { font-size: 10px; font-weight: 600; } .legend { font-size: 12px; }</style>',
    ]
    parts.append(f'<text x="{left}" y="34" class="title">Latent Map: Countries in Place, Arrows to Displaced Positions</text>')
    parts.append(f'<text x="{left}" y="52" class="subtitle">Model: {_svg_escape(model_label)}. Fixed anchors: {_svg_escape(", ".join(anchor_codes))}. Country surfaces stay in their real locations.</text>')
    for record in records:
        path_d = _rings_to_path(record["rings"], project)  # type: ignore[arg-type]
        parts.append(
            f'<path d="{path_d}" fill="{_fill_color(record, shifted=False)}" stroke="{_stroke_color(record, shifted=False)}" stroke-width="0.55" fill-rule="evenodd" />'
        )
    for record in records:
        start_x, start_y = project(float(record["centroid_lon"]), float(record["centroid_lat"]))
        end_x, end_y = project(float(record["target_lon"]), float(record["target_lat"]))
        color = _arrow_color(record)
        width_px = 1.6 if bool(record["is_spain_relevant"]) or bool(record["is_anchor"]) else 0.8
        opacity = 0.88 if bool(record["is_spain_relevant"]) or bool(record["is_anchor"]) else 0.34
        parts.append(f'<line x1="{start_x:.2f}" y1="{start_y:.2f}" x2="{end_x:.2f}" y2="{end_y:.2f}" stroke="{color}" stroke-width="{width_px}" opacity="{opacity}" marker-end="url(#arrowhead)" />')
        parts.append(f'<circle cx="{start_x:.2f}" cy="{start_y:.2f}" r="1.8" fill="#46525b" opacity="0.55" />')
        parts.append(f'<circle cx="{end_x:.2f}" cy="{end_y:.2f}" r="2.4" fill="{color}" opacity="0.92" />')
        if bool(record["is_spain_relevant"]) or bool(record["is_anchor"]):
            parts.append(f'<text x="{end_x + 4:.2f}" y="{end_y - 4:.2f}" class="label">{_svg_escape(record['iso3'])}</text>')
    _draw_legend(parts, left, top, shifted=False)
    parts.append('</svg>')
    text = "\n".join(parts).replace("'</f'text>", "</text>")
    path.write_text(text, encoding="utf-8")


def _write_shifted_shapes_map(path: Path, records: list[dict[str, object]], anchor_codes: tuple[str, ...], model_label: str) -> None:
    width = 1960
    height = 1080
    left = 40
    top = 60
    bounds = _compute_bounds(records)
    project = _project_factory(bounds, width, height, left, top)
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#f8f3eb" />',
        '<defs><marker id="arrowhead" markerWidth="8" markerHeight="8" refX="7" refY="3.5" orient="auto"><polygon points="0 0, 8 3.5, 0 7" fill="#8c4545" /></marker></defs>',
        '<style>text { font-family: Segoe UI, Arial, sans-serif; fill: #243039; } .title { font-size: 28px; font-weight: 700; } .subtitle { font-size: 14px; fill: #57636c; } .label { font-size: 10px; font-weight: 600; } .legend { font-size: 12px; }</style>',
    ]
    parts.append(f'<text x="{left}" y="34" class="title">Latent Map: Shifted Country Shapes</text>')
    parts.append(f'<text x="{left}" y="52" class="subtitle">Model: {_svg_escape(model_label)}. Fixed anchors: {_svg_escape(", ".join(anchor_codes))}. Shapes are translated by each country\'s latent displacement vector.</text>')
    for record in records:
        path_d = _rings_to_path(record["rings"], project)  # type: ignore[arg-type]
        parts.append(f'<path d="{path_d}" fill="#f3f1ea" stroke="#d3d5d0" stroke-width="0.45" opacity="0.55" fill-rule="evenodd" />')
    for record in records:
        shifted_path_d = _rings_to_path(
            record["rings"], project, delta_lon=float(record["delta_lon"]), delta_lat=float(record["delta_lat"])
        )  # type: ignore[arg-type]
        parts.append(
            f'<path d="{shifted_path_d}" fill="{_fill_color(record, shifted=True)}" stroke="{_stroke_color(record, shifted=True)}" stroke-width="0.6" opacity="0.92" fill-rule="evenodd" />'
        )
    for record in records:
        start_x, start_y = project(float(record["centroid_lon"]), float(record["centroid_lat"]))
        end_x, end_y = project(float(record["target_lon"]), float(record["target_lat"]))
        color = _arrow_color(record)
        width_px = 1.6 if bool(record["is_spain_relevant"]) or bool(record["is_anchor"]) else 0.8
        opacity = 0.9 if bool(record["is_spain_relevant"]) or bool(record["is_anchor"]) else 0.28
        parts.append(f'<line x1="{start_x:.2f}" y1="{start_y:.2f}" x2="{end_x:.2f}" y2="{end_y:.2f}" stroke="{color}" stroke-width="{width_px}" opacity="{opacity}" marker-end="url(#arrowhead)" />')
        if bool(record["is_spain_relevant"]) or bool(record["is_anchor"]):
            parts.append(f'<text x="{end_x + 4:.2f}" y="{end_y - 4:.2f}" class="label">{_svg_escape(record['iso3'])}</text>')
    _draw_legend(parts, left, top, shifted=True)
    parts.append('</svg>')
    text = "\n".join(parts).replace("'</f'text>", "</text>")
    path.write_text(text, encoding="utf-8")


def build_latent_country_maps(settings: Settings, model_prefix: str = "latent_geodesic_ppml", output_dir: Path | None = None) -> list[Path]:
    records, anchor_codes, curated_base = _load_map_records(settings, model_prefix)
    output_dir = output_dir or settings.processed_dir / "models"
    output_dir.mkdir(parents=True, exist_ok=True)
    static_map_path = output_dir / f"{model_prefix}_country_arrows_map.svg"
    shifted_map_path = output_dir / f"{model_prefix}_country_shifted_shapes_map.svg"
    _write_static_arrows_map(static_map_path, records, anchor_codes, model_prefix)
    _write_shifted_shapes_map(shifted_map_path, records, anchor_codes, model_prefix)

    curated_dir = settings.root_dir / "output"
    curated_dir.mkdir(parents=True, exist_ok=True)
    curated_static_path = curated_dir / f"{curated_base}_country_arrows_map.svg"
    curated_shifted_path = curated_dir / f"{curated_base}_country_shifted_shapes_map.svg"
    shutil.copyfile(static_map_path, curated_static_path)
    shutil.copyfile(shifted_map_path, curated_shifted_path)
    return [static_map_path, shifted_map_path, curated_static_path, curated_shifted_path]
