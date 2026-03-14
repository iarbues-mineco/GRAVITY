from __future__ import annotations

import re
from pathlib import Path

from .harmonize import harmonize_source_frame
from .settings import Settings


WORKBOOK_PATTERNS = [
    "destination",
    "origin",
]


def _require_pandas():
    try:
        import pandas as pd  # type: ignore
    except ImportError as exc:  # pragma: no cover
        raise RuntimeError("pandas is required for UN DESA stock normalization commands.") from exc
    return pd


def _require_openpyxl() -> None:
    try:
        import openpyxl  # noqa: F401
    except ImportError as exc:  # pragma: no cover
        raise RuntimeError(
            "openpyxl is required to read the UN DESA workbook. Install project dependencies again with `python -m pip install -e .`."
        ) from exc


def _normalize_text(value: object) -> str:
    return re.sub(r"\s+", " ", str(value or "").strip().lower())


def _resolve_workbook_path(settings: Settings, raw_path: Path | None = None) -> Path:
    if raw_path is not None:
        if not raw_path.exists():
            raise FileNotFoundError(f"UN DESA stock workbook not found: {raw_path}")
        return raw_path

    source_dir = settings.raw_dir / "un_desa_stock"
    if not source_dir.exists():
        raise FileNotFoundError(
            "UN DESA stock folder is missing. Place the `Destination and origin` workbook in `data/raw/un_desa_stock/` first."
        )

    workbooks = sorted(source_dir.glob("*.xlsx"))
    if not workbooks:
        raise FileNotFoundError(
            "No UN DESA stock workbook found. Place the `Destination and origin` workbook in `data/raw/un_desa_stock/` first."
        )

    preferred = [
        path
        for path in workbooks
        if all(token in path.stem.lower() for token in WORKBOOK_PATTERNS)
    ]
    if len(preferred) == 1:
        return preferred[0]
    if len(preferred) > 1:
        preferred_names = ", ".join(path.name for path in preferred)
        raise ValueError(
            "Multiple UN DESA workbooks look like destination-origin matrices. "
            f"Please keep only one or pass an explicit path. Candidates: {preferred_names}"
        )
    if len(workbooks) == 1:
        return workbooks[0]

    workbook_names = ", ".join(path.name for path in workbooks)
    raise ValueError(
        "Multiple UN DESA workbooks found, but none was uniquely identifiable as the destination-origin matrix. "
        f"Files present: {workbook_names}"
    )


def _find_destination_origin_sheet(raw_path: Path):
    pd = _require_pandas()

    excel = pd.ExcelFile(raw_path, engine="openpyxl")
    required_markers = ["location code of destination", "location code of origin"]

    for sheet_name in excel.sheet_names:
        preview = pd.read_excel(raw_path, sheet_name=sheet_name, header=None, nrows=120, engine="openpyxl")
        for row_index in range(len(preview.index)):
            row_text = " | ".join(_normalize_text(value) for value in preview.iloc[row_index].tolist())
            if all(marker in row_text for marker in required_markers):
                table = pd.read_excel(raw_path, sheet_name=sheet_name, header=row_index, engine="openpyxl")
                table = table.dropna(axis=0, how="all").dropna(axis=1, how="all")
                return table

    sheet_names = ", ".join(excel.sheet_names)
    raise ValueError(
        "Could not find a destination-origin data sheet in the UN DESA workbook. "
        f"Sheets inspected: {sheet_names}"
    )


def _find_column(columns, required_terms: list[str], forbidden_terms: list[str] | None = None, optional: bool = False):
    forbidden_terms = forbidden_terms or []
    normalized = {column: _normalize_text(column) for column in columns}
    for column, text in normalized.items():
        if all(term in text for term in required_terms) and not any(term in text for term in forbidden_terms):
            return column
    if optional:
        return None
    required_label = ", ".join(required_terms)
    raise ValueError(f"Could not find a UN DESA column containing: {required_label}")


def _parse_year_columns(columns) -> list[tuple[object, int, str]]:
    year_columns: list[tuple[object, int, str]] = []
    pattern = re.compile(r"^(?P<year>\d{4})(?:\.(?P<suffix>\d+))?$")

    for column in columns:
        text = str(column).strip()
        match = pattern.match(text)
        if not match:
            continue
        year = int(match.group("year"))
        suffix = match.group("suffix")
        if suffix is None:
            sex = "both_sexes"
        elif suffix == "1":
            sex = "male"
        elif suffix == "2":
            sex = "female"
        else:
            continue
        year_columns.append((column, year, sex))

    if not year_columns:
        raise ValueError("Could not identify any stock year columns in the UN DESA workbook.")
    return year_columns


def _coerce_stock_series(series):
    pd = _require_pandas()

    cleaned = (
        series.astype("string")
        .str.strip()
        .replace(
            {
                "": pd.NA,
                ".": pd.NA,
                "..": pd.NA,
                "...": pd.NA,
                "\u2026": pd.NA,
                "~": "0",
                "-": "0",
                "\u2013": "0",
                "\u2014": "0",
            }
        )
        .str.replace(",", "", regex=False)
        .str.replace(" ", "", regex=False)
    )
    return pd.to_numeric(cleaned, errors="coerce")


def _reshape_destination_origin_table(table):
    pd = _require_pandas()

    destination_name_column = _find_column(
        table.columns,
        required_terms=["destination"],
        forbidden_terms=["location code", "notes", "type of data"],
    )
    destination_code_column = _find_column(table.columns, required_terms=["location code of destination"])
    destination_notes_column = _find_column(table.columns, required_terms=["notes of destination"], optional=True)
    destination_type_column = _find_column(table.columns, required_terms=["type of data of destination"], optional=True)
    origin_name_column = _find_column(
        table.columns,
        required_terms=["origin"],
        forbidden_terms=["location code", "notes", "type of data"],
    )
    origin_code_column = _find_column(table.columns, required_terms=["location code of origin"])
    origin_type_column = _find_column(table.columns, required_terms=["type of data of origin"], optional=True)

    year_columns = _parse_year_columns(table.columns)

    metadata_columns = {
        "destination_name_raw": destination_name_column,
        "destination_code_raw": destination_code_column,
        "origin_name_raw": origin_name_column,
        "origin_code_raw": origin_code_column,
    }
    if destination_notes_column is not None:
        metadata_columns["destination_notes"] = destination_notes_column
    if destination_type_column is not None:
        metadata_columns["destination_data_type"] = destination_type_column
    if origin_type_column is not None:
        metadata_columns["origin_data_type"] = origin_type_column

    metadata = table[list(metadata_columns.values())].copy().rename(columns={value: key for key, value in metadata_columns.items()})
    metadata["destination_code_raw"] = metadata["destination_code_raw"].astype("string").str.strip()
    metadata["origin_code_raw"] = metadata["origin_code_raw"].astype("string").str.strip()
    metadata["destination_name_raw"] = (
        metadata["destination_name_raw"].astype("string").str.replace(r"\*+$", "", regex=True).str.strip()
    )
    metadata["origin_name_raw"] = (
        metadata["origin_name_raw"].astype("string").str.replace(r"\*+$", "", regex=True).str.strip()
    )

    long_frames = []
    for source_column, stock_year, sex in year_columns:
        subset = metadata.copy()
        subset["stock_year"] = stock_year
        subset["sex"] = sex
        subset["migrant_stock"] = _coerce_stock_series(table[source_column])
        long_frames.append(subset)

    long_frame = pd.concat(long_frames, ignore_index=True)
    long_frame = long_frame.loc[
        long_frame[["destination_name_raw", "origin_name_raw", "destination_code_raw", "origin_code_raw"]].notna().any(axis=1)
    ].copy()

    wide = (
        long_frame.pivot_table(
            index=[column for column in long_frame.columns if column not in {"sex", "migrant_stock"}],
            columns="sex",
            values="migrant_stock",
            aggfunc="first",
        )
        .reset_index()
    )
    wide.columns.name = None

    for sex in ["both_sexes", "male", "female"]:
        if sex not in wide.columns:
            wide[sex] = pd.NA

    wide = wide.rename(
        columns={
            "both_sexes": "migrant_stock_both_sexes",
            "male": "migrant_stock_male",
            "female": "migrant_stock_female",
        }
    )
    wide = wide.loc[
        wide[["migrant_stock_both_sexes", "migrant_stock_male", "migrant_stock_female"]].notna().any(axis=1)
    ].copy()
    wide["source_dataset"] = "un_desa_stock"
    wide["stock_unit"] = "persons"
    wide["stock_reference"] = "mid_year"
    return wide


def _build_entity_map(frame, settings: Settings, code_column: str, name_column: str, prefix: str):
    unique_entities = frame[[code_column, name_column]].drop_duplicates().copy()
    unique_entities[code_column] = unique_entities[code_column].astype("string")
    unique_entities[name_column] = unique_entities[name_column].astype("string")

    harmonized, unmatched = harmonize_source_frame(
        unique_entities,
        settings,
        source_id="un_desa_stock",
        code_column=code_column,
        name_column=name_column,
    )

    mapping = harmonized[
        [
            code_column,
            name_column,
            "canonical_id",
            "iso3",
            "canonical_name",
            "entity_type",
            "is_current",
            "match_method",
        ]
    ].copy()
    mapping = mapping.rename(
        columns={
            "canonical_id": f"{prefix}_canonical_id",
            "iso3": f"{prefix}_iso3",
            "canonical_name": f"{prefix}_name",
            "entity_type": f"{prefix}_entity_type",
            "is_current": f"{prefix}_is_current",
            "match_method": f"{prefix}_match_method",
        }
    )

    counts = frame.groupby([code_column, name_column], dropna=False).size().reset_index(name="row_count")
    unmatched = unmatched[[code_column, name_column]].drop_duplicates().merge(counts, on=[code_column, name_column], how="left")
    unmatched = unmatched.rename(
        columns={
            code_column: f"{prefix}_code_raw",
            name_column: f"{prefix}_name_raw",
        }
    )
    return mapping, unmatched


def normalize_un_desa_stock(
    settings: Settings,
    raw_path: Path | None = None,
    output_dir: Path | None = None,
) -> list[Path]:
    _require_openpyxl()
    pd = _require_pandas()

    workbook_path = _resolve_workbook_path(settings, raw_path=raw_path)
    table = _find_destination_origin_sheet(workbook_path)
    stock = _reshape_destination_origin_table(table)

    destination_map, destination_unmatched = _build_entity_map(
        stock,
        settings,
        code_column="destination_code_raw",
        name_column="destination_name_raw",
        prefix="destination",
    )
    origin_map, origin_unmatched = _build_entity_map(
        stock,
        settings,
        code_column="origin_code_raw",
        name_column="origin_name_raw",
        prefix="origin",
    )

    normalized = stock.merge(destination_map, on=["destination_code_raw", "destination_name_raw"], how="left")
    normalized = normalized.merge(origin_map, on=["origin_code_raw", "origin_name_raw"], how="left")
    normalized["is_self_corridor"] = normalized["origin_canonical_id"].eq(normalized["destination_canonical_id"])

    for column in ["migrant_stock_both_sexes", "migrant_stock_male", "migrant_stock_female"]:
        normalized[column] = pd.to_numeric(normalized[column], errors="coerce")

    output_dir = output_dir or settings.processed_dir / "stocks"
    output_dir.mkdir(parents=True, exist_ok=True)

    stock_path = output_dir / "bilateral_migrant_stock_un_desa.csv"
    origin_unmatched_path = output_dir / "un_desa_stock_unmatched_origins.csv"
    destination_unmatched_path = output_dir / "un_desa_stock_unmatched_destinations.csv"

    output_columns = [
        "source_dataset",
        "stock_year",
        "stock_reference",
        "stock_unit",
        "destination_code_raw",
        "destination_name_raw",
        "destination_notes",
        "destination_data_type",
        "destination_canonical_id",
        "destination_iso3",
        "destination_name",
        "destination_entity_type",
        "destination_is_current",
        "destination_match_method",
        "origin_code_raw",
        "origin_name_raw",
        "origin_data_type",
        "origin_canonical_id",
        "origin_iso3",
        "origin_name",
        "origin_entity_type",
        "origin_is_current",
        "origin_match_method",
        "migrant_stock_both_sexes",
        "migrant_stock_male",
        "migrant_stock_female",
        "is_self_corridor",
    ]
    existing_output_columns = [column for column in output_columns if column in normalized.columns]
    normalized = normalized[existing_output_columns].sort_values(
        ["stock_year", "origin_code_raw", "destination_code_raw"]
    ).reset_index(drop=True)

    normalized.to_csv(stock_path, index=False)
    origin_unmatched.to_csv(origin_unmatched_path, index=False)
    destination_unmatched.to_csv(destination_unmatched_path, index=False)
    return [stock_path, origin_unmatched_path, destination_unmatched_path]
