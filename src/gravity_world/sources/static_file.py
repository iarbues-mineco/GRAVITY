from __future__ import annotations

from pathlib import Path

from ..catalog import SourceSpec
from ..http import download_to_file


class StaticFileSource:
    def __init__(self, spec: SourceSpec, raw_dir: Path) -> None:
        self.spec = spec
        self.base_dir = raw_dir / spec.output_subdir
        self.base_dir.mkdir(parents=True, exist_ok=True)

    def download(self) -> Path:
        if not self.spec.url or not self.spec.filename:
            raise ValueError(f"Static source {self.spec.id} is missing url or filename.")
        return download_to_file(self.spec.url, self.base_dir / self.spec.filename)
