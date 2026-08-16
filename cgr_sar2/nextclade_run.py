"""Nextclade CLI: dataset atualizado + atribuicao Pango / clado Nextstrain."""

from __future__ import annotations

import csv
import io
import os
import platform
import stat
import subprocess
import tempfile
from dataclasses import dataclass, field
from pathlib import Path
from urllib.request import Request, urlopen

from .who_current import who_status_for_pango

DATASET_NAME = "nextstrain/sars-cov-2/wuhan-hu-1/orfs"
DEFAULT_DATASET_DIR = Path("data/nextclade/sars-cov-2")
DEFAULT_BIN_DIR = Path("bin")

_ASSET = {
    ("Linux", "x86_64"): "nextclade-x86_64-unknown-linux-gnu",
    ("Linux", "aarch64"): "nextclade-aarch64-unknown-linux-gnu",
    ("Darwin", "x86_64"): "nextclade-x86_64-apple-darwin",
    ("Darwin", "arm64"): "nextclade-aarch64-apple-darwin",
}


@dataclass
class PangoAssignment:
    header: str
    pango: str
    clade: str
    who_status: str
    who_label: str
    qc: str
    substitutions: int | None
    missing: int | None
    aa_substitutions: str
    extra: dict[str, str] = field(default_factory=dict)
    dataset_dir: str = ""
    error: str | None = None


def nextclade_binary(bin_dir: Path = DEFAULT_BIN_DIR) -> Path:
    explicit = os.environ.get("NEXTCLADE_BIN")
    if explicit:
        path = Path(explicit)
        if path.is_file():
            return path
    for name in ("nextclade", "nextclade.exe"):
        candidate = bin_dir / name
        if candidate.is_file():
            return candidate
    from shutil import which

    found = which("nextclade")
    if found:
        return Path(found)
    raise FileNotFoundError(
        "Nextclade nao encontrado. Execute `python -m cgr_sar2 update-nextclade` "
        "ou coloque o binario em bin/nextclade."
    )


def _asset_name() -> str:
    system = platform.system()
    machine = platform.machine()
    key = (system, machine)
    if key not in _ASSET:
        raise RuntimeError(f"Sistema nao suportado para Nextclade: {system} {machine}")
    return _ASSET[key]


def install_nextclade(bin_dir: Path = DEFAULT_BIN_DIR) -> Path:
    bin_dir.mkdir(parents=True, exist_ok=True)
    dest = bin_dir / "nextclade"
    url = (
        "https://github.com/nextstrain/nextclade/releases/latest/download/"
        + _asset_name()
    )
    req = Request(url, headers={"User-Agent": "cgr-sar2"})
    with urlopen(req, timeout=120) as resp:
        dest.write_bytes(resp.read())
    dest.chmod(dest.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
    return dest


def update_dataset(
    dataset_dir: Path = DEFAULT_DATASET_DIR,
    bin_dir: Path = DEFAULT_BIN_DIR,
    name: str = DATASET_NAME,
) -> Path:
    try:
        binary = nextclade_binary(bin_dir)
    except FileNotFoundError:
        binary = install_nextclade(bin_dir)
    dataset_dir.mkdir(parents=True, exist_ok=True)
    cmd = [
        str(binary),
        "dataset",
        "get",
        "--name",
        name,
        "--output-dir",
        str(dataset_dir),
    ]
    subprocess.run(cmd, check=True, capture_output=True, text=True)
    return dataset_dir


def ensure_ready(
    dataset_dir: Path = DEFAULT_DATASET_DIR,
    bin_dir: Path = DEFAULT_BIN_DIR,
) -> tuple[Path, Path]:
    try:
        binary = nextclade_binary(bin_dir)
    except FileNotFoundError:
        binary = install_nextclade(bin_dir)
    marker = dataset_dir / "pathogen.json"
    if not marker.is_file():
        update_dataset(dataset_dir=dataset_dir, bin_dir=bin_dir)
    return binary, dataset_dir


def assign_fasta(
    fasta_text: str,
    dataset_dir: Path = DEFAULT_DATASET_DIR,
    bin_dir: Path = DEFAULT_BIN_DIR,
) -> list[PangoAssignment]:
    binary, dataset = ensure_ready(dataset_dir=dataset_dir, bin_dir=bin_dir)
    with tempfile.TemporaryDirectory(prefix="cgr-sar2-nc-") as tmp:
        fasta_path = Path(tmp) / "query.fasta"
        tsv_path = Path(tmp) / "nextclade.tsv"
        fasta_path.write_text(fasta_text, encoding="utf-8")
        cmd = [
            str(binary),
            "run",
            "--input-dataset",
            str(dataset),
            "--output-tsv",
            str(tsv_path),
            "--silent",
            str(fasta_path),
        ]
        proc = subprocess.run(cmd, capture_output=True, text=True)
        if proc.returncode != 0 or not tsv_path.is_file():
            detail = (proc.stderr or proc.stdout or "Nextclade falhou").strip()
            raise RuntimeError(detail[:2000])
        text = tsv_path.read_text(encoding="utf-8", errors="replace")
    return _parse_tsv(text, dataset_dir=str(dataset))


def _parse_tsv(text: str, dataset_dir: str) -> list[PangoAssignment]:
    reader = csv.DictReader(io.StringIO(text), delimiter="\t")
    rows: list[PangoAssignment] = []
    for raw in reader:
        pango = (
            raw.get("Nextclade_pango")
            or raw.get("pangoLineage")
            or raw.get("NextcladePango")
            or ""
        ).strip()
        clade = (raw.get("clade") or raw.get("clade_nextstrain") or "").strip()
        status, label = who_status_for_pango(pango, clade)
        rows.append(
            PangoAssignment(
                header=raw.get("seqName") or raw.get("seqName") or "",
                pango=pango or "unassigned",
                clade=clade or "—",
                who_status=status,
                who_label=label,
                qc=(raw.get("qc.overallStatus") or raw.get("qcStatus") or "").strip(),
                substitutions=_maybe_int(raw.get("totalSubstitutions")),
                missing=_maybe_int(raw.get("totalMissing")),
                aa_substitutions=raw.get("aaSubstitutions") or raw.get("aaSubstitutions") or "",
                extra={
                    "partiallyAliased": raw.get("partiallyAliased") or "",
                    "qcScore": raw.get("qc.overallScore") or "",
                    "frameShifts": raw.get("frameShifts") or "",
                },
                dataset_dir=dataset_dir,
            )
        )
    if not rows:
        raise RuntimeError("Nextclade nao devolveu resultados. Verifique o FASTA.")
    return rows


def _maybe_int(value: str | None) -> int | None:
    if value is None or value == "":
        return None
    try:
        return int(float(value))
    except ValueError:
        return None


def dataset_stamp(dataset_dir: Path = DEFAULT_DATASET_DIR) -> str:
    pathogen = dataset_dir / "pathogen.json"
    if not pathogen.is_file():
        return "dataset nao instalado"
    try:
        import json

        data = json.loads(pathogen.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return str(pathogen.stat().st_mtime)
    meta = data.get("meta") or data.get("attributes") or {}
    version = (
        data.get("version")
        or meta.get("updatedAt")
        or meta.get("version")
        or data.get("schemaVersion")
    )
    if isinstance(version, dict):
        version = version.get("tag") or version.get("updatedAt")
    return str(version or "ok")
