"""Frequencias recentes de linhagens Pango via CovSpectrum LAPIS (GenBank aberto)."""

from __future__ import annotations

import json
from datetime import date, timedelta
from urllib.error import URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen

LAPIS_AGG = "https://lapis.cov-spectrum.org/open/v2/sample/aggregated"
DEFAULT_DAYS = 60
MAX_LINEAGES = 15


def recent_pango_frequencies(
    days: int = DEFAULT_DAYS,
    limit: int = MAX_LINEAGES,
) -> dict:
    """Agrega amostras publicas dos ultimos `days` dias por pangoLineage."""
    end = date.today()
    start = end - timedelta(days=max(days, 14))
    params = {
        "fields": "pangoLineage",
        "dateFrom": start.isoformat(),
        "dateTo": end.isoformat(),
        "limit": 500,
    }
    url = f"{LAPIS_AGG}?{urlencode(params)}"
    req = Request(url, headers={"Accept": "application/json", "User-Agent": "cgr-sar2"})
    try:
        with urlopen(req, timeout=30) as resp:
            payload = json.loads(resp.read().decode("utf-8"))
    except (URLError, TimeoutError, json.JSONDecodeError, OSError) as exc:
        return {
            "ok": False,
            "error": str(exc),
            "start": start.isoformat(),
            "end": end.isoformat(),
            "rows": [],
            "total": 0,
        }
    rows = payload.get("data") or []
    cleaned = []
    for row in rows:
        lineage = row.get("pangoLineage") or "unassigned"
        count = int(row.get("count") or 0)
        if count <= 0:
            continue
        cleaned.append({"pango": lineage, "count": count})
    cleaned.sort(key=lambda item: item["count"], reverse=True)
    total = sum(item["count"] for item in cleaned)
    top = cleaned[:limit]
    for item in top:
        item["percent"] = 100.0 * item["count"] / total if total else 0.0
    return {
        "ok": True,
        "start": start.isoformat(),
        "end": end.isoformat(),
        "source": "CovSpectrum LAPIS open (GenBank)",
        "rows": top,
        "total": total,
    }
