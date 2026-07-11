"""Load and represent project configuration from ``config.yaml``."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import yaml


def _project_root() -> Path:
    """Repository root (two levels up from this file: src/mirna_tcga/)."""
    return Path(__file__).resolve().parents[2]


@dataclass
class Config:
    """Typed view over the YAML config with convenient id-builder helpers."""

    cbioportal: dict[str, Any]
    studies: dict[str, str]
    profiles: dict[str, str]
    sample_lists: dict[str, str]
    analysis: dict[str, Any]
    cohorts: dict[str, list[str]] = field(default_factory=dict)
    xena: dict[str, Any] = field(default_factory=dict)
    xena_mirna: dict[str, str] = field(default_factory=dict)
    raw: dict[str, Any] = field(default_factory=dict, repr=False)

    # -- id builders -------------------------------------------------------
    def mrna_profile(self, study_key: str) -> str:
        return self.studies[study_key] + self.profiles["mrna_suffix"]

    def mutation_profile(self, study_key: str) -> str:
        return self.studies[study_key] + self.profiles["mutation_suffix"]

    def cna_profile(self, study_key: str) -> str:
        return self.studies[study_key] + self.profiles["cna_suffix"]

    def mirna_profile(self, study_key: str) -> str:
        return self.studies[study_key] + self.profiles["mirna_suffix"]

    def all_samples_list(self, study_key: str) -> str:
        return self.studies[study_key] + self.sample_lists["all_suffix"]

    def cna_samples_list(self, study_key: str) -> str:
        return self.studies[study_key] + self.sample_lists["cna_suffix"]

    def sequenced_samples_list(self, study_key: str) -> str:
        return self.studies[study_key] + self.sample_lists["sequenced_suffix"]


def load_config(path: str | Path | None = None) -> Config:
    """Read ``config.yaml`` (defaults to the one at the repo root)."""
    path = Path(path) if path is not None else _project_root() / "config.yaml"
    with open(path, "r", encoding="utf-8") as fh:
        data = yaml.safe_load(fh)
    return Config(
        cbioportal=data["cbioportal"],
        studies=data["studies"],
        profiles=data["profiles"],
        sample_lists=data["sample_lists"],
        analysis=data["analysis"],
        cohorts=data.get("cohorts", {}),
        xena=data.get("xena", {}),
        xena_mirna=data.get("xena_mirna", {}),
        raw=data,
    )
