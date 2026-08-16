"""Gradio interface for variant identification."""

from __future__ import annotations

from functools import lru_cache
from pathlib import Path

from .labels import VARIANT_LABELS, VARIANT_WHO_STATUS
from .predict import VariantClassifier, fcgr_preview_rgb

PAPER = (
    "Camara, G. B. M., Coutinho, M. G. F., Barbosa, R. de M., & Fernandes, M. A. C. "
    "(2026). Integrating Chaos Game Representation with Convolutional Neural Networks "
    "for Accurate Classification of SARS-CoV-2 Variants. IEEE Access. "
    "DOI: 10.1109/ACCESS.2026.3723929"
)


@lru_cache(maxsize=1)
def _classifier(weights: str) -> VariantClassifier:
    return VariantClassifier(weights_path=weights or None)


def classify(fasta_file, sequence_text, weights_path: str | None = None):
    path = None
    if fasta_file is not None:
        path = fasta_file if isinstance(fasta_file, (str, Path)) else fasta_file.name
    clf = _classifier(str(weights_path or ""))
    preds = clf.predict_fasta(fasta_text=sequence_text or None, fasta_path=path)
    pred = preds[0]
    extra = f"\n\nOther sequences in the file: {len(preds) - 1}" if len(preds) > 1 else ""
    warn = f"\n\n**Warning:** {pred.warning}" if pred.warning else ""
    summary = (
        f"### {pred.variant}\n"
        f"{VARIANT_WHO_STATUS[pred.variant]}  \n"
        f"Confidence: **{pred.confidence:.2%}**  \n"
        f"Length: {pred.length:,} nt · ambiguous bases: {pred.ambiguous_fraction:.2%}"
        f"{warn}{extra}"
    )
    probs = pred.probabilities
    labels = list(VARIANT_LABELS)
    image = fcgr_preview_rgb(pred.fcgr)
    table = [[name, VARIANT_WHO_STATUS[name], f"{probs[name]:.4f}"] for name in labels]
    return summary, image, table


def launch(
    weights_path: str | Path | None = None,
    server_name: str = "127.0.0.1",
    server_port: int = 7860,
    share: bool = False,
) -> None:
    import gradio as gr

    weights = str(weights_path) if weights_path else ""

    def _run(fasta_file, sequence_text):
        try:
            return classify(fasta_file, sequence_text, weights)
        except FileNotFoundError as exc:
            raise gr.Error(str(exc)) from exc
        except ValueError as exc:
            raise gr.Error(str(exc)) from exc

    with gr.Blocks(title="CGR-SAR2 - SARS-CoV-2 variants") as demo:
        gr.Markdown(
            """
# CGR-SAR2
Identificacao de variantes de **SARS-CoV-2** a partir do genoma completo,
usando Frequency Chaos Game Representation (FCGR, k = 6) e a CNN do artigo.

O modelo distingue **Alpha, Beta, Gamma, Delta, Iota e Epsilon**.
Sequencias de outras linhagens (por exemplo Omicron) serao forcadas a uma destas classes.
Interprete a confianca com cuidado.
"""
        )
        with gr.Row():
            with gr.Column():
                fasta = gr.File(label="Ficheiro FASTA", file_types=[".fasta", ".fa", ".fna", ".txt"])
                sequence = gr.Textbox(
                    label="Ou cole a sequencia / FASTA",
                    lines=8,
                    placeholder=">EPI_ISL_...\nATGTAG...",
                )
                btn = gr.Button("Identificar variante", variant="primary")
            with gr.Column():
                summary = gr.Markdown()
                image = gr.Image(label="Imagem FCGR 64x64", type="numpy")
                table = gr.Dataframe(headers=["Variante", "WHO", "Probabilidade"], interactive=False)
        btn.click(_run, inputs=[fasta, sequence], outputs=[summary, image, table])
        gr.Markdown(f"*Referencia:* {PAPER}")

    demo.launch(server_name=server_name, server_port=server_port, share=share)


if __name__ == "__main__":
    launch()
