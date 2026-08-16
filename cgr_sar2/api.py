"""API HTTP: POST /classify."""

from __future__ import annotations

from typing import Literal

from fastapi import FastAPI, File, Form, HTTPException, UploadFile
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel, Field

from .classify import classify_fasta_text, engine_status
from .nextclade_run import dataset_stamp, ensure_ready

app = FastAPI(
    title="CGR-SAR2",
    description="Atribuicao Pango (Nextclade e, opcionalmente, Pangolin) para SARS-CoV-2.",
    version="1.1.0",
)
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_methods=["*"],
    allow_headers=["*"],
)


class ClassifyBody(BaseModel):
    fasta: str = Field(..., min_length=10)
    engine: Literal["nextclade", "pangolin", "both"] = "nextclade"


@app.on_event("startup")
def _startup() -> None:
    ensure_ready()


@app.get("/health")
def health() -> dict:
    return {"ok": True, "dataset": dataset_stamp(), "engines": engine_status()}


def _run(text: str, engine: str) -> dict:
    if engine not in {"nextclade", "pangolin", "both"}:
        raise HTTPException(400, "engine deve ser nextclade, pangolin ou both")
    try:
        rows = classify_fasta_text(text, engine=engine)
    except FileNotFoundError as exc:
        raise HTTPException(503, str(exc)) from exc
    except Exception as exc:
        raise HTTPException(422, str(exc)[:2000]) from exc
    return {"n": len(rows), "engine": engine, "results": rows}


@app.post("/classify")
def classify_json(body: ClassifyBody) -> dict:
    return _run(body.fasta, body.engine)


@app.post("/classify/file")
async def classify_file(
    fasta: UploadFile = File(...),
    engine: str = Form("nextclade"),
) -> dict:
    text = (await fasta.read()).decode("utf-8", errors="replace")
    if not text.strip():
        raise HTTPException(400, "Ficheiro FASTA vazio.")
    return _run(text, engine)
