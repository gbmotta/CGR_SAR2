"""Download a NCBI Datasets virus genome zip and keep N quality FASTA files."""

from __future__ import annotations

import urllib.parse
import urllib.request
import zipfile
from pathlib import Path

from .cgr import ambiguous_fraction
from .fasta import parse_fasta
from .labels import MAX_AMBIGUOUS_FRACTION, MIN_GENOME_LENGTH

DOWNLOAD_URL = "https://api.ncbi.nlm.nih.gov/datasets/v2/virus/taxon/2697049/genome/download"


def _quality_ok(sequence: str) -> bool:
    return len(sequence) >= MIN_GENOME_LENGTH and ambiguous_fraction(sequence) <= MAX_AMBIGUOUS_FRACTION


def download_lineage_package(
    lineage: str,
    dest_dir: Path,
    n_keep: int,
    extra_params: dict | None = None,
    zip_path: Path | None = None,
) -> int:
    dest_dir.mkdir(parents=True, exist_ok=True)
    existing = list(dest_dir.glob("*.fasta"))
    if len(existing) >= n_keep:
        print(f"{lineage}: {len(existing)} files already present")
        return len(existing)

    params = {
        "pangolin_classification": lineage,
        "complete_only": "true",
        "filename": f"{lineage}.zip",
    }
    if extra_params:
        params.update(extra_params)
    url = DOWNLOAD_URL + "?" + urllib.parse.urlencode(params)
    zip_path = zip_path or Path("data/tmp") / f"{lineage.replace('.', '_')}.zip"
    zip_path.parent.mkdir(parents=True, exist_ok=True)
    print(f"{lineage}: downloading package to {zip_path} ...")
    req = urllib.request.Request(url, headers={"User-Agent": "cgr-sar2", "Accept": "application/zip"})
    with urllib.request.urlopen(req, timeout=600) as response, zip_path.open("wb") as out:
        while True:
            chunk = response.read(1024 * 1024)
            if not chunk:
                break
            out.write(chunk)
    print(f"{lineage}: zip size {zip_path.stat().st_size / 1e6:.1f} MB")

    seen = {p.stem for p in existing}
    kept = len(existing)
    with zipfile.ZipFile(zip_path) as zf:
        fasta_names = [n for n in zf.namelist() if n.endswith((".fna", ".fasta", ".fa"))]
        print(f"{lineage}: fasta members {fasta_names[:8]}")
        for name in fasta_names:
            if kept >= n_keep:
                break
            with zf.open(name) as handle:
                text = handle.read().decode("utf-8", errors="replace")
            for header, seq in parse_fasta(text):
                if kept >= n_keep:
                    break
                seq_u = seq.upper()
                if not _quality_ok(seq_u):
                    continue
                acc = header.split()[0].replace("/", "_")
                if acc in seen:
                    continue
                seen.add(acc)
                (dest_dir / f"{acc}.fasta").write_text(f">{header}\n{seq_u}\n", encoding="utf-8")
                kept += 1
            print(f"  {lineage}: {kept}/{n_keep}")
    try:
        zip_path.unlink()
    except OSError:
        pass
    return kept
