"""Tabela de resultados para CSV/XLSX e API."""

from __future__ import annotations

import csv
import io
from pathlib import Path

from .nextclade_run import PangoAssignment
from .signatures import best_signature, match_signatures, signature_summary


def assignment_to_row(
    asn: PangoAssignment,
    length: int | None = None,
    pangolin: str | None = None,
    pangolin_qc: str | None = None,
) -> dict:
    best, score = best_signature(asn.aa_substitutions)
    sigs = match_signatures(asn.aa_substitutions)
    row = {
        "header": asn.header,
        "pango_nextclade": asn.pango,
        "pango_pangolin": pangolin or "",
        "agreement": (
            bool(pangolin)
            and pangolin.split()[0] == asn.pango.split()[0]
            if pangolin
            else ""
        ),
        "nextstrain": asn.clade,
        "who": asn.who_status,
        "who_label": asn.who_label,
        "qc_nextclade": asn.qc,
        "qc_pangolin": pangolin_qc or "",
        "length": length if length is not None else "",
        "substitutions": asn.substitutions if asn.substitutions is not None else "",
        "missing": asn.missing if asn.missing is not None else "",
        "signature_best": best,
        "signature_score": round(score, 3),
        "signatures": signature_summary(asn.aa_substitutions),
        "aa_substitutions": asn.aa_substitutions,
    }
    for name, info in sigs.items():
        row[f"sig_{name}_present"] = ",".join(info["present"])
        row[f"sig_{name}_missing"] = ",".join(info["missing"])
    return row


def rows_to_csv(rows: list[dict]) -> str:
    if not rows:
        return ""
    buf = io.StringIO()
    writer = csv.DictWriter(buf, fieldnames=list(rows[0].keys()))
    writer.writeheader()
    writer.writerows(rows)
    return buf.getvalue()


def write_csv(rows: list[dict], path: Path) -> Path:
    path.write_text(rows_to_csv(rows), encoding="utf-8")
    return path


def write_xlsx(rows: list[dict], path: Path) -> Path:
    import pandas as pd

    pd.DataFrame(rows).to_excel(path, index=False)
    return path
