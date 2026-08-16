"""Pangolin / UShER opcional (binario local ou imagem Docker staphb/pangolin)."""

from __future__ import annotations

import csv
import os
import shutil
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path


DEFAULT_IMAGE = os.environ.get("PANGOLIN_DOCKER", "staphb/pangolin:latest")


@dataclass
class PangolinAssignment:
    header: str
    pango: str
    qc: str
    scorpio: str
    version: str
    note: str
    error: str | None = None


def pangolin_available() -> tuple[bool, str]:
    if shutil.which("pangolin"):
        return True, "binary"
    if shutil.which("docker") and os.environ.get("CGR_SAR2_PANGOLIN", "1") != "0":
        return True, "docker"
    return False, "unavailable"


def assign_pangolin(fasta_text: str) -> list[PangolinAssignment]:
    mode_ok, mode = pangolin_available()
    if not mode_ok:
        raise FileNotFoundError(
            "Pangolin nao encontrado. Instale o pacote ou use Docker: "
            "docker pull staphb/pangolin:latest"
        )
    with tempfile.TemporaryDirectory(prefix="cgr-sar2-pg-") as tmp:
        work = Path(tmp)
        fasta_path = work / "query.fasta"
        out_path = work / "lineage_report.csv"
        fasta_path.write_text(fasta_text, encoding="utf-8")
        if mode == "binary":
            cmd = ["pangolin", str(fasta_path), "--outfile", str(out_path)]
            proc = subprocess.run(cmd, capture_output=True, text=True)
        else:
            cmd = [
                "docker",
                "run",
                "--rm",
                "-v",
                f"{work}:/data",
                DEFAULT_IMAGE,
                "pangolin",
                "/data/query.fasta",
                "--outfile",
                "/data/lineage_report.csv",
            ]
            proc = subprocess.run(cmd, capture_output=True, text=True)
        if proc.returncode != 0 or not out_path.is_file():
            detail = (proc.stderr or proc.stdout or "Pangolin falhou").strip()
            raise RuntimeError(detail[:2000])
        return _parse_report(out_path.read_text(encoding="utf-8", errors="replace"))


def _parse_report(text: str) -> list[PangolinAssignment]:
    reader = csv.DictReader(io_string(text))
    rows: list[PangolinAssignment] = []
    for raw in reader:
        rows.append(
            PangolinAssignment(
                header=(raw.get("taxon") or raw.get("Sequence name") or "").strip(),
                pango=(raw.get("lineage") or raw.get("Lineage") or "").strip(),
                qc=(raw.get("qc_status") or raw.get("status") or "").strip(),
                scorpio=(raw.get("scorpio_call") or "").strip(),
                version=(raw.get("pangolin_version") or raw.get("version") or "").strip(),
                note=(raw.get("note") or raw.get("qc_notes") or "").strip(),
            )
        )
    if not rows:
        raise RuntimeError("Pangolin nao devolveu resultados.")
    return rows


def io_string(text: str):
    import io

    return io.StringIO(text)
