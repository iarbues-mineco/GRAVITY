from __future__ import annotations

from pathlib import Path

from ..catalog import IndicatorSpec, SourceSpec
from ..http import build_url, fetch_json
from ..io import write_csv, write_json


class WorldBankSource:
    def __init__(self, spec: SourceSpec, raw_dir: Path) -> None:
        self.spec = spec
        self.base_dir = raw_dir / spec.output_subdir
        self.base_dir.mkdir(parents=True, exist_ok=True)

    def download(self) -> list[Path]:
        outputs = [self._download_country_metadata()]
        for indicator in self.spec.indicators:
            outputs.append(self._download_indicator(indicator))
        return outputs

    def _download_country_metadata(self) -> Path:
        all_rows: list[dict] = []
        page = 1
        total_pages = 1

        while page <= total_pages:
            url = build_url(
                self.spec.country_metadata_url or "",
                {
                    "format": self.spec.format or "json",
                    "per_page": self.spec.per_page or 500,
                    "page": page,
                },
            )
            payload = fetch_json(url)
            meta, rows = payload
            total_pages = int(meta["pages"])
            all_rows.extend(rows)
            page += 1

        json_path = self.base_dir / "country_metadata.json"
        write_json(json_path, all_rows)

        csv_rows = []
        for row in all_rows:
            csv_rows.append(
                {
                    "id": row.get("id"),
                    "iso2c": row.get("iso2Code"),
                    "name": row.get("name"),
                    "region_id": (row.get("region") or {}).get("id"),
                    "region_name": (row.get("region") or {}).get("value"),
                    "income_level_id": (row.get("incomeLevel") or {}).get("id"),
                    "income_level_name": (row.get("incomeLevel") or {}).get("value"),
                    "lending_type_id": (row.get("lendingType") or {}).get("id"),
                    "lending_type_name": (row.get("lendingType") or {}).get("value"),
                    "capital_city": row.get("capitalCity"),
                    "longitude": row.get("longitude"),
                    "latitude": row.get("latitude"),
                }
            )

        return write_csv(
            self.base_dir / "country_metadata.csv",
            csv_rows,
            fieldnames=[
                "id",
                "iso2c",
                "name",
                "region_id",
                "region_name",
                "income_level_id",
                "income_level_name",
                "lending_type_id",
                "lending_type_name",
                "capital_city",
                "longitude",
                "latitude",
            ],
        )

    def _download_indicator(self, indicator: IndicatorSpec) -> Path:
        all_rows: list[dict] = []
        page = 1
        total_pages = 1

        while page <= total_pages:
            url = build_url(
                (self.spec.base_url or "").format(indicator=indicator.code),
                {
                    "format": self.spec.format or "json",
                    "per_page": self.spec.per_page or 20000,
                    "page": page,
                    "date": f"{indicator.date_from}:{indicator.date_to}",
                },
            )
            payload = fetch_json(url)
            meta, rows = payload
            total_pages = int(meta["pages"])
            all_rows.extend(rows)
            page += 1

        normalized_rows = []
        for row in all_rows:
            normalized_rows.append(
                {
                    "indicator_code": (row.get("indicator") or {}).get("id"),
                    "indicator_name": (row.get("indicator") or {}).get("value"),
                    "country_name": (row.get("country") or {}).get("value"),
                    "country_iso3": row.get("countryiso3code"),
                    "year": row.get("date"),
                    "value": row.get("value"),
                    "unit": row.get("unit"),
                    "obs_status": row.get("obs_status"),
                    "decimal": row.get("decimal"),
                }
            )

        return write_csv(
            self.base_dir / f"{indicator.slug}.csv",
            normalized_rows,
            fieldnames=[
                "indicator_code",
                "indicator_name",
                "country_name",
                "country_iso3",
                "year",
                "value",
                "unit",
                "obs_status",
                "decimal",
            ],
        )
