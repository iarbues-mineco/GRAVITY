from __future__ import annotations

import json
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import Request, urlopen


USER_AGENT = "gravity-world-migration/0.1"


def build_url(base_url: str, params: dict[str, str | int]) -> str:
    query = urlencode(params)
    separator = "&" if "?" in base_url else "?"
    return f"{base_url}{separator}{query}"


def fetch_json(url: str, timeout: int = 60) -> list | dict:
    request = Request(url, headers={"User-Agent": USER_AGENT})
    with urlopen(request, timeout=timeout) as response:
        return json.loads(response.read().decode("utf-8"))


def download_to_file(url: str, target_path: Path, timeout: int = 120) -> Path:
    target_path.parent.mkdir(parents=True, exist_ok=True)
    request = Request(url, headers={"User-Agent": USER_AGENT})
    with urlopen(request, timeout=timeout) as response, target_path.open("wb") as output_file:
        while True:
            chunk = response.read(1024 * 1024)
            if not chunk:
                break
            output_file.write(chunk)
    return target_path
