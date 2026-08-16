"""Estatuto WHO e nomes amigaveis para linhagens Pango atuais (julho 2026)."""

from __future__ import annotations

# Fonte: https://www.who.int/activities/tracking-SARS-CoV-2-variants
# VUMs em circulacao (27 Jul 2026) e VOI JN.1. Prefixos mais especificos primeiro.
_WHO_PREFIXES: tuple[tuple[str, str, str], ...] = (
    ("PQ.16.1.1", "VUM", "PQ.16.1.1 (descendente de NB.1.8.1)"),
    ("PQ.", "VUM", "PQ.* (familia NB.1.8.1 / 25B)"),
    ("NB.1.8.1", "VUM", "NB.1.8.1 (clado Nextstrain 25B)"),
    ("NB.", "VUM", "NB.* (familia 25B)"),
    ("XFG", "VUM", "XFG (clado Nextstrain 25C)"),
    ("BA.3.2", "VUM", "BA.3.2"),
    ("KP.", "VOI", "JN.1 / Omicron (alias KP.*)"),
    ("LP.", "VOI", "JN.1 / Omicron (alias LP.*)"),
    ("KS.", "VOI", "JN.1 / Omicron (alias KS.*)"),
    ("XEC", "VOI", "XEC (recombinante da era JN.1)"),
    ("JN.1", "VOI", "JN.1 / Omicron (clado Nextstrain 24A)"),
    ("BA.2.86", "VOI", "BA.2.86 / Omicron (ancestral de JN.1)"),
    ("XBB", "historica", "Omicron XBB"),
    ("BA.5", "historica", "Omicron BA.5"),
    ("BA.4", "historica", "Omicron BA.4"),
    ("BA.2", "historica", "Omicron BA.2"),
    ("BA.1", "historica", "Omicron BA.1"),
    ("B.1.1.529", "historica", "Omicron"),
    ("AY.", "historica", "Delta (AY.*)"),
    ("B.1.617.2", "historica", "VOC Delta"),
    ("B.1.1.7", "historica", "VOC Alpha"),
    ("Q.", "historica", "Alpha (Q.*)"),
    ("B.1.351", "historica", "VOC Beta"),
    ("P.1", "historica", "VOC Gamma"),
    ("B.1.526", "historica", "VOI Iota"),
    ("B.1.427", "historica", "VOI Epsilon"),
    ("B.1.429", "historica", "VOI Epsilon"),
)


def _matches(lineage: str, prefix: str) -> bool:
    lin = lineage.strip()
    if not lin or lin in {".", "unassigned", "Unassigned"}:
        return False
    if prefix.endswith("."):
        return lin == prefix[:-1] or lin.startswith(prefix)
    return lin == prefix or lin.startswith(prefix + ".")


def who_status_for_pango(lineage: str | None, clade: str | None = None) -> tuple[str, str]:
    """Devolve (estatuto, descricao) para uma linhagem Pango e clado Nextstrain."""
    if not lineage:
        return _from_clade(clade, "Linhagem Pango nao atribuida")
    text = lineage.strip()
    if text.lower() in {"", ".", "unassigned", "none", "nan"}:
        return _from_clade(clade, "Linhagem Pango nao atribuida")
    for prefix, status, label in _WHO_PREFIXES:
        if _matches(text, prefix):
            return status, f"{label} · Pango {text}"
    clade_hit = _from_clade(clade, None)
    if clade_hit is not None:
        status, label = clade_hit
        return status, f"{label} · Pango {text}"
    if text.startswith("BA.") or text.startswith("BQ.") or text.startswith("CH."):
        return "historica", f"Omicron / descendente · Pango {text}"
    return "sem designacao WHO", f"Linhagem Pango {text} (sem letra grega / VUM propria)"


def _from_clade(clade: str | None, fallback: str | None) -> tuple[str, str] | None:
    if not clade or clade in {"—", "-", "."}:
        if fallback is None:
            return None
        return "indeterminado", fallback
    tag = clade.strip().upper()
    if tag.startswith("25C"):
        return "VUM", f"Clado Nextstrain {clade} (era XFG)"
    if tag.startswith("25B"):
        return "VUM", f"Clado Nextstrain {clade} (era NB.1.8.1)"
    if tag.startswith("24"):
        return "VOI", f"Clado Nextstrain {clade} (era JN.1)"
    if fallback is None:
        return None
    return "indeterminado", fallback


def circulating_highlights() -> list[dict[str, str]]:
    """Resumo estatico das variantes WHO em circulacao (julho 2026)."""
    return [
        {
            "pango": "JN.1",
            "nextstrain": "24A",
            "who": "VOI",
            "nota": "VOI em circulacao; muitas linhagens atuais sao descendentes (KP.*, LP.*, NB.*, XFG, ...)",
        },
        {
            "pango": "XFG",
            "nextstrain": "25C",
            "who": "VUM",
            "nota": "Recombinante JN.1; uma das mais frequentes em 2026",
        },
        {
            "pango": "NB.1.8.1",
            "nextstrain": "25B",
            "who": "VUM",
            "nota": "Predominante em algumas regioes (ex. WPRO)",
        },
        {
            "pango": "PQ.16.1.1",
            "nextstrain": "25B",
            "who": "VUM",
            "nota": "Designada em 27 Jul 2026 (NB.1.8.1 + mutacoes S)",
        },
        {
            "pango": "BA.3.2",
            "nextstrain": "—",
            "who": "VUM",
            "nota": "Linha BA.3 distinta, monitorizada desde Dez 2025",
        },
    ]
