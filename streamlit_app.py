"""App Streamlit: lote Pango, assinaturas WHO e exportacao."""

from __future__ import annotations

import streamlit as st

from cgr_sar2.cgr import ambiguous_fraction, sequence_to_fcgr
from cgr_sar2.circulation import recent_pango_frequencies
from cgr_sar2.classify import classify_records, engine_status
from cgr_sar2.fasta import records_from_input
from cgr_sar2.nextclade_run import dataset_stamp, ensure_ready
from cgr_sar2.predict import fcgr_preview_rgb
from cgr_sar2.report import rows_to_csv, write_xlsx
from cgr_sar2.signatures import SPIKE_SIGNATURES, match_signatures
from cgr_sar2.who_current import circulating_highlights

st.set_page_config(page_title="CGR-SAR2 · Pango atual", page_icon="🧬", layout="wide")


@st.cache_resource(show_spinner="A preparar Nextclade e o dataset SARS-CoV-2…")
def _prepare_nextclade() -> str:
    ensure_ready()
    return dataset_stamp()


@st.cache_data(ttl=6 * 3600, show_spinner=False)
def _circulation(days: int = 60):
    return recent_pango_frequencies(days=days)


def _records_from_uploads(files) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    if not files:
        return records
    for uploaded in files:
        text = uploaded.getvalue().decode("utf-8", errors="replace")
        records.extend(records_from_input(text, None))
    return records


def main() -> None:
    st.title("CGR-SAR2 · identificação atualizada")
    st.caption(
        "Linhagem **Pango** (Nextclade; Pangolin/UShER opcional em Docker) e clado **Nextstrain**. "
        "O CGR+CNN do artigo (Alpha–Epsilon) é só histórico."
    )

    try:
        stamp = _prepare_nextclade()
        engines = engine_status()
        st.success(f"Dataset Nextclade `{stamp}` · Pangolin: {engines['pangolin']}")
    except Exception as exc:
        st.error(
            "Nextclade indisponível. Local: `python -m cgr_sar2 update-nextclade`. "
            "Na Cloud, confirme o `setup.sh`."
        )
        st.exception(exc)
        return

    left, right = st.columns((1.2, 1.0), gap="large")

    with left:
        st.subheader("Genomas")
        uploaded = st.file_uploader(
            "Um ou mais FASTA",
            type=["fasta", "fa", "fna", "txt"],
            accept_multiple_files=True,
        )
        pasted = st.text_area(
            "Ou cole FASTA / nucleótidos (várias sequências)",
            height=160,
            placeholder=">amostra_1\nATTAAA...\n>amostra_2\nATTAAA...",
        )
        engine_choice = "nextclade"
        if engines["pangolin"] != "unavailable":
            engine_choice = st.radio(
                "Motor",
                ["nextclade", "both"],
                horizontal=True,
                help="both corre Nextclade e Pangolin (UShER) e compara as linhagens.",
            )
        run = st.button("Identificar linhagens", type="primary", use_container_width=True)

        if run:
            records = _records_from_uploads(uploaded)
            records.extend(records_from_input(pasted or None, None))
            if not records:
                st.warning("Forneça FASTA ou cole sequências.")
            else:
                try:
                    rows = classify_records(records, engine=engine_choice)
                except Exception as exc:
                    st.error("Falha na classificação.")
                    st.exception(exc)
                    return
                st.session_state["rows"] = rows
                st.session_state["records"] = records

        rows = st.session_state.get("rows")
        records = st.session_state.get("records") or []
        if rows:
            st.dataframe(rows, hide_index=True, use_container_width=True)
            csv_text = rows_to_csv(rows)
            st.download_button("Descarregar CSV", csv_text, "cgr_sar2_pango.csv", "text/csv")
            try:
                from pathlib import Path
                import tempfile

                tmp = Path(tempfile.gettempdir()) / "cgr_sar2_pango.xlsx"
                write_xlsx(rows, tmp)
                st.download_button(
                    "Descarregar Excel",
                    tmp.read_bytes(),
                    "cgr_sar2_pango.xlsx",
                    "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                )
            except Exception:
                pass
            for rec, row in zip(records, rows):
                _render_detail(rec[0], rec[1], row)

    with right:
        st.subheader("O que circula agora")
        st.markdown("WHO (julho 2026) e frequências públicas GenBank (LAPIS, 60 dias).")
        st.dataframe(circulating_highlights(), hide_index=True, use_container_width=True)
        circ = _circulation(60)
        if circ.get("ok") and circ.get("rows"):
            st.caption(f"{circ['total']:,} genomas · {circ['start']} a {circ['end']}")
            st.bar_chart(
                {"Pango": [r["pango"] for r in circ["rows"]], "%": [round(r["percent"], 2) for r in circ["rows"]]},
                x="Pango",
                y="%",
                horizontal=True,
            )
        st.markdown("**Assinaturas da espícula (WHO)**")
        for name, muts in SPIKE_SIGNATURES.items():
            st.caption(f"{name}: {', '.join(muts)}")


def _render_detail(header: str, sequence: str, row: dict) -> None:
    st.markdown(f"#### {header or '(sem cabeçalho)'}")
    c1, c2, c3, c4 = st.columns(4)
    c1.metric("Pango Nextclade", row.get("pango_nextclade") or "—")
    c2.metric("Nextstrain", row.get("nextstrain") or "—")
    c3.metric("WHO", row.get("who") or "—")
    c4.metric("Assinatura", f"{row.get('signature_best')} ({row.get('signature_score')})")
    if row.get("pango_pangolin"):
        st.write(f"Pangolin/UShER: **{row['pango_pangolin']}** · acordo: {row.get('agreement')}")
    st.write(row.get("who_label") or "")
    st.markdown(
        f"Comprimento **{len(sequence):,} nt** · N **{ambiguous_fraction(sequence):.2%}** · "
        f"QC {row.get('qc_nextclade') or '—'}"
    )
    sigs = match_signatures(row.get("aa_substitutions") or "")
    sig_table = [
        {
            "variante": name,
            "presentes": ", ".join(info["present"]) or "—",
            "em falta": ", ".join(info["missing"]) or "—",
            "cobertura": f"{info['n_present']}/{info['n_total']}",
        }
        for name, info in sigs.items()
    ]
    st.dataframe(sig_table, hide_index=True, use_container_width=True)
    if row.get("aa_substitutions"):
        with st.expander("Todas as substituições de aminoácidos"):
            st.text(row["aa_substitutions"])
    try:
        st.image(
            fcgr_preview_rgb(sequence_to_fcgr(sequence)),
            caption="CGR 64×64 (visualização)",
        )
    except Exception:
        pass
    st.divider()


main()
