"""Classificacao unificada Nextclade e, se disponivel, Pangolin."""

from __future__ import annotations

from .fasta import parse_fasta, records_from_input, records_from_paths
from .nextclade_run import assign_fasta
from .pangolin_run import PangolinAssignment, assign_pangolin, pangolin_available
from .report import assignment_to_row


def collect_records(fasta_text: str | None, paths: list | None) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    if paths:
        records.extend(records_from_paths(paths))
    if fasta_text and fasta_text.strip():
        records.extend(records_from_input(fasta_text, None))
    return records


def classify_records(
    records: list[tuple[str, str]],
    engine: str = "nextclade",
) -> list[dict]:
    if not records:
        raise ValueError("Nenhuma sequencia FASTA.")
    blob = "\n".join(f">{header}\n{seq}" for header, seq in records)
    lengths = {header.split()[0]: len(seq) for header, seq in records}
    lengths.update({header: len(seq) for header, seq in records})

    nc_rows = assign_fasta(blob)
    pg_map: dict[str, PangolinAssignment] = {}
    if engine in {"pangolin", "both"}:
        for item in assign_pangolin(blob):
            pg_map[_key(item.header)] = item
    if engine == "pangolin" and not pg_map:
        raise RuntimeError("Pangolin nao devolveu linhagens.")

    out: list[dict] = []
    for asn in nc_rows:
        key = _key(asn.header)
        pg = pg_map.get(key)
        length = lengths.get(asn.header) or lengths.get(key)
        out.append(
            assignment_to_row(
                asn,
                length=length,
                pangolin=pg.pango if pg else None,
                pangolin_qc=pg.qc if pg else None,
            )
        )
    return out


def classify_fasta_text(fasta_text: str, engine: str = "nextclade") -> list[dict]:
    records = list(parse_fasta(fasta_text) if fasta_text.strip().startswith(">") else [])
    if not records and fasta_text.strip():
        records = [("query", "".join(fasta_text.split()))]
    return classify_records(records, engine=engine)


def _key(header: str) -> str:
    return (header or "").split()[0]


def engine_status() -> dict[str, str]:
    ok, mode = pangolin_available()
    return {
        "nextclade": "ready",
        "pangolin": mode if ok else "unavailable",
    }
