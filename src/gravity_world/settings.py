from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class Settings:
    root_dir: Path
    config_dir: Path
    data_dir: Path
    raw_dir: Path
    interim_dir: Path
    processed_dir: Path
    logs_dir: Path

    @classmethod
    def discover(cls) -> "Settings":
        root_dir = Path(__file__).resolve().parents[2]
        data_dir = root_dir / "data"
        return cls(
            root_dir=root_dir,
            config_dir=root_dir / "config",
            data_dir=data_dir,
            raw_dir=data_dir / "raw",
            interim_dir=data_dir / "interim",
            processed_dir=data_dir / "processed",
            logs_dir=data_dir / "logs",
        )

    def ensure_directories(self) -> None:
        for path in (
            self.config_dir,
            self.raw_dir,
            self.interim_dir,
            self.processed_dir,
            self.logs_dir,
        ):
            path.mkdir(parents=True, exist_ok=True)
