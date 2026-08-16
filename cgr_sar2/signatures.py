"""Mutacoes-assinatura da espícula (WHO VUM / VOI, julho 2026)."""

from __future__ import annotations

# Fonte: https://www.who.int/activities/tracking-SARS-CoV-2-variants
SPIKE_SIGNATURES: dict[str, tuple[str, ...]] = {
    "XFG": (
        "S:T22N",
        "S:S31P",
        "S:K182R",
        "S:R190S",
        "S:R346T",
        "S:K444R",
        "S:V445R",
        "S:F456L",
        "S:N487D",
        "S:Q493E",
        "S:T572I",
    ),
    "NB.1.8.1": (
        "S:T22N",
        "S:F59S",
        "S:G184S",
        "S:A435S",
        "S:F456L",
        "S:T478I",
        "S:Q493E",
    ),
    "PQ.16.1.1": (
        "S:D253G",
        "S:N417T",
        "S:D420N",
        "S:I478T",
    ),
    "JN.1": ("S:L455S",),
}


def parse_aa_mutations(aa_substitutions: str) -> set[str]:
    found: set[str] = set()
    if not aa_substitutions:
        return found
    for raw in aa_substitutions.replace(";", ",").split(","):
        token = raw.strip()
        if not token:
            continue
        if token.startswith("Spike:"):
            token = "S:" + token.split(":", 1)[1]
        found.add(token)
    return found


def match_signatures(aa_substitutions: str) -> dict[str, dict]:
    present = parse_aa_mutations(aa_substitutions)
    report: dict[str, dict] = {}
    for name, muts in SPIKE_SIGNATURES.items():
        hit = [m for m in muts if m in present]
        miss = [m for m in muts if m not in present]
        report[name] = {
            "present": hit,
            "missing": miss,
            "score": len(hit) / len(muts) if muts else 0.0,
            "n_present": len(hit),
            "n_total": len(muts),
        }
    return report


def signature_summary(aa_substitutions: str) -> str:
    parts = []
    for name, row in match_signatures(aa_substitutions).items():
        parts.append(f"{name} {row['n_present']}/{row['n_total']}")
    return "; ".join(parts)


def best_signature(aa_substitutions: str) -> tuple[str, float]:
    ranked = match_signatures(aa_substitutions)
    if not ranked:
        return "—", 0.0
    covered = [key for key, row in ranked.items() if row["score"] >= 0.6]
    pool = covered or list(ranked)
    name = max(pool, key=lambda key: (ranked[key]["n_present"], ranked[key]["score"]))
    return name, float(ranked[name]["score"])
