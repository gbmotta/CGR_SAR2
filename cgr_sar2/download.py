"""Download public SARS-CoV-2 genomes from NCBI Datasets + E-utilities."""

from __future__ import annotations

import argparse
import http.client
import io
import json
import random
import time
import urllib.error
import urllib.parse
import urllib.request
import zipfile
from pathlib import Path

from .cgr import ambiguous_fraction
from .fasta import parse_fasta
from .labels import (
    MAX_AMBIGUOUS_FRACTION,
    MIN_GENOME_LENGTH,
    NCBI_LINEAGES,
    VARIANT_LABELS,
)

REPORT_URL = "https://api.ncbi.nlm.nih.gov/datasets/v2/virus/taxon/2697049/genome/download"
EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
DEFAULT_OUT = Path("data/fasta")
FETCH_BATCH = 20
HUMAN_TAXID = 9606


def _http_get(url: str, retries: int = 6, timeout: int = 180) -> bytes:
    last_error: Exception | None = None
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "cgr-sar2/1.0"})
            with urllib.request.urlopen(req, timeout=timeout) as response:
                return response.read()
        except (urllib.error.URLError, TimeoutError, http.client.IncompleteRead) as exc:
            last_error = exc
            time.sleep(min(30, 2 ** attempt))
    raise RuntimeError(f"HTTP failure for {url}: {last_error}")


def _is_human(record: dict) -> bool:
    host = record.get("host") or {}
    lineage = host.get("lineage") or []
    for item in lineage:
        if str(item.get("taxId") or item.get("tax_id") or "") == str(HUMAN_TAXID):
            return True
        if str(item.get("name") or "").lower() == "homo sapiens":
            return True
    return False


def list_accessions(lineage: str, extra_params: dict | None = None) -> list[str]:
    params = {
        "pangolin_classification": lineage,
        "complete_only": "true",
        "filename": f"{lineage}.zip",
    }
    if extra_params:
        params.update(extra_params)
    url = REPORT_URL + "?" + urllib.parse.urlencode(params)
    print(f"  fetching accession report for {lineage}...")
    blob = _http_get(url, timeout=300)
    accessions: list[str] = []
    with zipfile.ZipFile(io.BytesIO(blob)) as zf:
        report_name = next((n for n in zf.namelist() if n.endswith("data_report.jsonl")), None)
        if report_name is None:
            print(f"  no data_report.jsonl in package for {lineage}")
            return []
        with zf.open(report_name) as handle:
            for raw in handle:
                rec = json.loads(raw)
                if rec.get("completeness") != "COMPLETE":
                    continue
                length = int(rec.get("length") or 0)
                if length < MIN_GENOME_LENGTH:
                    continue
                if rec.get("host") and not _is_human(rec):
                    continue
                acc = rec.get("accession")
                if acc:
                    accessions.append(str(acc))
    print(f"  {lineage}: {len(accessions)} complete human genomes in report")
    return accessions


def fetch_fasta(accessions: list[str]) -> str:
    chunks: list[str] = []
    for start in range(0, len(accessions), FETCH_BATCH):
        batch = accessions[start : start + FETCH_BATCH]
        params = {
            "db": "nuccore",
            "id": ",".join(batch),
            "rettype": "fasta",
            "retmode": "text",
        }
        url = EFETCH_URL + "?" + urllib.parse.urlencode(params)
        chunks.append(_http_get(url).decode("utf-8", errors="replace"))
        time.sleep(0.4)
    return "\n".join(chunks)


def _quality_ok(sequence: str) -> bool:
    return len(sequence) >= MIN_GENOME_LENGTH and ambiguous_fraction(sequence) <= MAX_AMBIGUOUS_FRACTION


def download_variant(variant: str, out_dir: Path, per_class: int, seed: int) -> int:
    dest = out_dir / variant.lower()
    dest.mkdir(parents=True, exist_ok=True)
    existing = list(dest.glob("*.fasta"))
    if len(existing) >= per_class:
        print(f"{variant}: {len(existing)} files already present, skipping.")
        return len(existing)

    rng = random.Random(seed)
    seen = {p.stem for p in existing}
    kept = len(existing)

    extra = {"usa_state": "CA"} if variant == "Epsilon" else None
    for lineage in NCBI_LINEAGES[variant]:
        if kept >= per_class:
            break
        print(f"{variant}: lineage {lineage}")
        try:
            ids = list_accessions(lineage, extra_params=extra)
        except RuntimeError as exc:
            print(f"  warning: report failed ({exc})")
            continue
        rng.shuffle(ids)
        for start in range(0, len(ids), FETCH_BATCH):
            if kept >= per_class:
                break
            batch = ids[start : start + FETCH_BATCH]
            try:
                fasta = fetch_fasta(batch)
            except RuntimeError as exc:
                print(f"  warning: batch failed ({exc})")
                continue
            if ">" not in fasta:
                print("  warning: NCBI did not return FASTA, skipping batch")
                continue
            for header, seq in parse_fasta(fasta):
                if kept >= per_class:
                    break
                seq_u = seq.upper()
                if not _quality_ok(seq_u):
                    continue
                acc = header.split()[0].replace("/", "_")
                if acc in seen:
                    continue
                seen.add(acc)
                (dest / f"{acc}.fasta").write_text(f">{header}\n{seq_u}\n", encoding="utf-8")
                kept += 1
            print(f"  {variant}: {kept}/{per_class} genomes saved")
    return kept


CURRENT_EXAMPLE_LINEAGES = (
    ("XFG", "xfg_example.fasta"),
    ("NB.1.8.1", "nb181_example.fasta"),
    ("JN.1", "jn1_example.fasta"),
)


def download_current_examples(out_dir: Path = Path("examples")) -> list[Path]:
    """Um genoma publico NCBI por linhagem atual (XFG, NB.1.8.1, JN.1)."""
    out_dir.mkdir(parents=True, exist_ok=True)
    written: list[Path] = []
    for lineage, filename in CURRENT_EXAMPLE_LINEAGES:
        dest = out_dir / filename
        print(f"Exemplo atual: {lineage}")
        try:
            ids = list_accessions(lineage)
        except RuntimeError as exc:
            print(f"  falhou o relatorio NCBI ({exc})")
            continue
        saved = False
        for start in range(0, min(len(ids), 40), FETCH_BATCH):
            batch = ids[start : start + FETCH_BATCH]
            try:
                fasta = fetch_fasta(batch)
            except RuntimeError as exc:
                print(f"  falhou o efetch ({exc})")
                continue
            for header, seq in parse_fasta(fasta):
                seq_u = seq.upper()
                if not _quality_ok(seq_u):
                    continue
                dest.write_text(f">{header}\n{seq_u}\n", encoding="utf-8")
                print(f"  gravado {dest} ({header.split()[0]}, {len(seq_u)} nt)")
                written.append(dest)
                saved = True
                break
            if saved:
                break
        if not saved:
            print(f"  nenhum genoma completo aceite para {lineage}")
    return written


def download_all(out_dir: Path, per_class: int, seed: int) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    for variant in VARIANT_LABELS:
        n = download_variant(variant, out_dir, per_class, seed)
        print(f"-> {variant}: {n} sequences")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Download NCBI genomes by WHO variant")
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--per-class", type=int, default=400)
    parser.add_argument("--seed", type=int, default=42)
    return parser


def main(argv: list[str] | None = None) -> None:
    args = build_parser().parse_args(argv)
    download_all(args.out_dir, args.per_class, args.seed)


if __name__ == "__main__":
    main()
