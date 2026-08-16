"""FASTA parsing."""

from __future__ import annotations

from collections.abc import Iterator
from pathlib import Path


def _looks_like_path(source: str) -> bool:
    if "\n" in source or source.startswith(">"):
        return False
    if len(source) > 4096:
        return False
    return Path(source).is_file()


def parse_fasta(source: str | Path) -> Iterator[tuple[str, str]]:
    if isinstance(source, Path) or (isinstance(source, str) and _looks_like_path(source)):
        text = Path(source).read_text(encoding="utf-8", errors="replace")
    else:
        text = str(source)
    header: str | None = None
    chunks: list[str] = []
    for raw_line in text.splitlines():
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                yield header, "".join(chunks)
            header = line[1:].strip()
            chunks = []
        else:
            chunks.append(line)
    if header is not None:
        yield header, "".join(chunks)


def records_from_input(fasta_text: str | None, fasta_path: str | Path | None) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    if fasta_path:
        records.extend(parse_fasta(Path(fasta_path)))
    if fasta_text and fasta_text.strip():
        stripped = fasta_text.strip()
        if stripped.startswith(">"):
            records.extend(parse_fasta(stripped))
        else:
            records.append(("pasted_sequence", "".join(stripped.split())))
    return records
