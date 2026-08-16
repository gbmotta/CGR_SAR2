"""App Streamlit: linhagem Pango atual (Nextclade) + circulacao recente."""

from __future__ import annotations

import streamlit as st

from cgr_sar2.cgr import ambiguous_fraction, sequence_to_fcgr
from cgr_sar2.circulation import recent_pango_frequencies
from cgr_sar2.fasta import records_from_input
from cgr_sar2.nextclade_run import assign_fasta, dataset_stamp, ensure_ready
from cgr_sar2.predict import fcgr_preview_rgb
from cgr_sar2.who_current import circulating_highlights

st.set_page_config(
    page_title="CGR-SAR2 · Pango atual",
    page_icon="🧬",
    layout="wide",
)


@st.cache_resource(show_spinner="A preparar Nextclade e o dataset SARS-CoV-2…")
def _prepare_nextclade() -> str:
    ensure_ready()
    return dataset_stamp()


@st.cache_data(ttl=6 * 3600, show_spinner=False)
def _circulation(days: int = 60):
    return recent_pango_frequencies(days=days)


def _fasta_from_upload(uploaded) -> str | None:
    if uploaded is None:
        return None
    return uploaded.getvalue().decode("utf-8", errors="replace")


def main() -> None:
    st.title("CGR-SAR2 · identificação atualizada")
    st.caption(
        "Atribuição **Pango** e clado **Nextstrain** com o dataset Nextclade mais recente. "
        "O classificador CGR+CNN do artigo (Alpha–Epsilon, 2021) fica só como contexto histórico."
    )

    try:
        stamp = _prepare_nextclade()
        st.success(f"Dataset Nextclade pronto · versão `{stamp}`")
    except Exception as exc:
        st.error(
            "Não foi possível instalar o Nextclade. No Streamlit Cloud, confirme que "
            "`setup.sh` correu. Localmente: `python -m cgr_sar2 update-nextclade`."
        )
        st.exception(exc)
        return

    left, right = st.columns((1.15, 1.0), gap="large")

    with left:
        st.subheader("Genoma")
        uploaded = st.file_uploader("FASTA", type=["fasta", "fa", "fna", "txt"])
        pasted = st.text_area(
            "Ou cole FASTA / nucleótidos",
            height=180,
            placeholder=">amostra\nATTAAAGGTTTATACCTTCCCAGGTAACAA...",
        )
        run = st.button("Identificar linhagem", type="primary", use_container_width=True)

        if run:
            fasta_text = _fasta_from_upload(uploaded)
            records = records_from_input(pasted or None, None)
            if fasta_text:
                records = records_from_input(fasta_text, None) + records
            if not records:
                st.warning("Forneça um ficheiro FASTA ou cole uma sequência.")
            else:
                blob = "\n".join(f">{h}\n{s}" for h, s in records)
                try:
                    assignments = assign_fasta(blob)
                except Exception as exc:
                    st.error("Nextclade não conseguiu analisar a sequência.")
                    st.exception(exc)
                    return
                for rec, asn in zip(records, assignments):
                    _render_result(rec[0], rec[1], asn)

    with right:
        st.subheader("O que circula agora")
        st.markdown(
            "Resumo WHO (julho 2026) e frequências **públicas** no GenBank "
            "(CovSpectrum LAPIS, últimos 60 dias)."
        )
        who_rows = circulating_highlights()
        st.dataframe(who_rows, hide_index=True, use_container_width=True)
        circ = _circulation(60)
        if circ.get("ok") and circ.get("rows"):
            st.caption(
                f"{circ['total']:,} genomas públicos de {circ['start']} a {circ['end']} · {circ['source']}"
            )
            chart = {
                "Pango": [r["pango"] for r in circ["rows"]],
                "%": [round(r["percent"], 2) for r in circ["rows"]],
            }
            st.bar_chart(chart, x="Pango", y="%", horizontal=True)
        elif circ.get("error"):
            st.info("Frequências em tempo real indisponíveis neste momento.")

        st.markdown(
            "Para **atualizar** o classificador: `python -m cgr_sar2 update-nextclade` "
            "ou faça *Reboot* da app no Streamlit Cloud (o `setup.sh` volta a descarregar o dataset)."
        )


def _render_result(header: str, sequence: str, asn) -> None:
    st.markdown(f"#### {header or '(sem cabeçalho)'}")
    c1, c2, c3, c4 = st.columns(4)
    c1.metric("Pango", asn.pango or "unassigned")
    c2.metric("Nextstrain", asn.clade)
    c3.metric("WHO", asn.who_status)
    c4.metric("QC Nextclade", asn.qc or "—")
    st.write(asn.who_label)
    meta = (
        f"Comprimento **{len(sequence):,} nt** · bases ambíguas "
        f"**{ambiguous_fraction(sequence):.2%}**"
    )
    if asn.substitutions is not None:
        meta += f" · substituições {asn.substitutions}"
    if asn.missing is not None:
        meta += f" · Ns {asn.missing}"
    st.markdown(meta)
    if asn.aa_substitutions:
        with st.expander("Substituições de aminoácidos (Nextclade)"):
            st.text(asn.aa_substitutions)
    try:
        image = fcgr_preview_rgb(sequence_to_fcgr(sequence))
        st.image(image, caption="CGR 64×64 (visualização; não é o classificador Pango)")
    except Exception:
        pass
    st.divider()


main()
