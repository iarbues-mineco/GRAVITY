from __future__ import annotations

from pathlib import Path

from ..catalog import SourceSpec
from ..http import download_to_file


class ManualFileSource:
    def __init__(self, spec: SourceSpec, raw_dir: Path) -> None:
        self.spec = spec
        self.base_dir = raw_dir / spec.output_subdir
        self.base_dir.mkdir(parents=True, exist_ok=True)

    def prepare(self) -> Path:
        if self.spec.url and self.spec.filename:
            return download_to_file(self.spec.url, self.base_dir / self.spec.filename)

        instructions = [
            f"Source id: {self.spec.id}",
            f"Description: {self.spec.description}",
            f"Landing page: {self.spec.landing_url or 'n/a'}",
            f"Expected filenames: {', '.join(self.spec.expected_filenames) if self.spec.expected_filenames else 'n/a'}",
            f"Notes: {self.spec.notes or 'n/a'}",
            "",
            "Action:",
            "1. Download the file manually from the landing page, or",
            "2. Add a direct `url` override in config/source_overrides.json and rerun the downloader.",
        ]
        instructions_path = self.base_dir / "MANUAL_DOWNLOAD.txt"
        instructions_path.write_text("\n".join(instructions), encoding="utf-8")
        return instructions_path
