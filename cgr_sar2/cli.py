"""Command-line interface: download, train, predict, app."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from .labels import VARIANT_LABELS
from .model import count_parameters


def _cmd_predict(args: argparse.Namespace) -> None:
    from .predict import VariantClassifier

    clf = VariantClassifier(weights_path=args.weights)
    results = clf.predict_fasta(fasta_text=args.sequence, fasta_path=args.fasta)
    payload = []
    for pred in results:
        row = {
            "header": pred.header,
            "variant": pred.variant,
            "who": pred.status,
            "confidence": round(pred.confidence, 6),
            "length": pred.length,
            "ambiguous_fraction": round(pred.ambiguous_fraction, 6),
            "probabilities": {k: round(v, 6) for k, v in pred.probabilities.items()},
            "warning": pred.warning,
        }
        payload.append(row)
        print(f"{pred.header or '(no header)'}")
        print(f"  Variant: {pred.variant} ({pred.status})")
        print(f"  Confidence: {pred.confidence:.2%}")
        if pred.warning:
            print(f"  Warning: {pred.warning}")
        print()
    if args.json:
        Path(args.json).write_text(json.dumps(payload, indent=2), encoding="utf-8")


def _cmd_train(args: argparse.Namespace) -> None:
    from .train import train_model

    train_model(
        data_dir=args.data_dir,
        cache_dir=args.cache_dir,
        out_path=args.out,
        epochs=args.epochs,
        batch_size=args.batch_size,
        learning_rate=args.lr,
        seed=args.seed,
        from_original=args.from_original,
        original_dir=args.original_dir,
        zscore=args.zscore,
    )


def _cmd_download(args: argparse.Namespace) -> None:
    from .download import download_all

    download_all(args.out_dir, args.per_class, args.seed)


def _cmd_app(args: argparse.Namespace) -> None:
    from .app import launch

    launch(weights_path=args.weights, server_name=args.host, server_port=args.port, share=args.share)


def _cmd_info(_args: argparse.Namespace) -> None:
    print(f"Classes historicas (artigo): {', '.join(VARIANT_LABELS)}")
    print(f"CNN parameters: {count_parameters():,}")
    print("Classificador atual: Nextclade + Pango (python -m cgr_sar2 pango)")


def _cmd_pango(args: argparse.Namespace) -> None:
    from .fasta import records_from_input
    from .nextclade_run import assign_fasta

    records = records_from_input(args.sequence, args.fasta)
    if not records:
        raise SystemExit("Forneca --fasta ou --sequence")
    blob = "\n".join(f">{header}\n{seq}" for header, seq in records)
    rows = assign_fasta(blob)
    payload = []
    for row in rows:
        item = {
            "header": row.header,
            "pango": row.pango,
            "nextstrain": row.clade,
            "who": row.who_status,
            "who_label": row.who_label,
            "qc": row.qc,
            "substitutions": row.substitutions,
            "missing": row.missing,
        }
        payload.append(item)
        print(row.header or "(sem cabecalho)")
        print(f"  Pango: {row.pango}")
        print(f"  Nextstrain: {row.clade}")
        print(f"  WHO: {row.who_status} - {row.who_label}")
        print(f"  QC: {row.qc or '-'}")
        print()
    if args.json:
        Path(args.json).write_text(json.dumps(payload, indent=2), encoding="utf-8")


def _cmd_update_nextclade(_args: argparse.Namespace) -> None:
    from .nextclade_run import dataset_stamp, update_dataset

    path = update_dataset()
    print(f"Dataset atualizado em {path}")
    print(f"Versao: {dataset_stamp(path)}")


def _cmd_streamlit(args: argparse.Namespace) -> None:
    import subprocess
    import sys

    app = Path(__file__).resolve().parent.parent / "streamlit_app.py"
    cmd = [
        sys.executable,
        "-m",
        "streamlit",
        "run",
        str(app),
        "--server.address",
        args.host,
        "--server.port",
        str(args.port),
    ]
    raise SystemExit(subprocess.call(cmd))


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="cgr-sar2",
        description="SARS-CoV-2: Pango/Nextclade atual + CNN historica do artigo.",
    )
    sub = parser.add_subparsers(dest="command", required=True)

    p_pred = sub.add_parser("predict", help="Classify a FASTA file or a pasted sequence")
    p_pred.add_argument("--fasta", type=Path, help="FASTA file (one or more sequences)")
    p_pred.add_argument("--sequence", help="Pasted sequence (FASTA or raw nucleotides)")
    p_pred.add_argument("--weights", type=Path, help="Path to variant_cnn.pt")
    p_pred.add_argument("--json", type=Path, help="Write results as JSON")
    p_pred.set_defaults(func=_cmd_predict)

    p_train = sub.add_parser("train", help="Train the CNN from original .mat images or FASTA")
    p_train.add_argument("--data-dir", type=Path, default=Path("data/fasta"))
    p_train.add_argument("--cache-dir", type=Path, default=Path("data/fcgr"))
    p_train.add_argument("--original-dir", type=Path, default=Path("data/original"))
    p_train.add_argument(
        "--from-original",
        action="store_true",
        help="Use GISAID CGR images from data/original (MATLAB Experimento 1)",
    )
    p_train.add_argument("--out", type=Path, default=Path("models/variant_cnn.pt"))
    p_train.add_argument("--epochs", type=int, default=40)
    p_train.add_argument("--batch-size", type=int, default=128)
    p_train.add_argument("--lr", type=float, default=1e-3)
    p_train.add_argument("--seed", type=int, default=42)
    p_train.add_argument("--zscore", action="store_true")
    p_train.set_defaults(func=_cmd_train)

    p_dl = sub.add_parser("download", help="Download public NCBI genomes")
    p_dl.add_argument("--out-dir", type=Path, default=Path("data/fasta"))
    p_dl.add_argument("--per-class", type=int, default=400)
    p_dl.add_argument("--seed", type=int, default=42)
    p_dl.set_defaults(func=_cmd_download)

    p_app = sub.add_parser("app", help="Open the Gradio web interface (modelo historico)")
    p_app.add_argument("--weights", type=Path)
    p_app.add_argument("--host", default="127.0.0.1")
    p_app.add_argument("--port", type=int, default=7860)
    p_app.add_argument("--share", action="store_true")
    p_app.set_defaults(func=_cmd_app)

    p_st = sub.add_parser("streamlit", help="Abrir a app Streamlit (Pango atual)")
    p_st.add_argument("--host", default="127.0.0.1")
    p_st.add_argument("--port", type=int, default=8501)
    p_st.set_defaults(func=_cmd_streamlit)

    p_pango = sub.add_parser("pango", help="Atribuir linhagem Pango com Nextclade")
    p_pango.add_argument("--fasta", type=Path)
    p_pango.add_argument("--sequence")
    p_pango.add_argument("--json", type=Path)
    p_pango.set_defaults(func=_cmd_pango)

    p_up = sub.add_parser("update-nextclade", help="Descarregar Nextclade e o dataset SARS-CoV-2")
    p_up.set_defaults(func=_cmd_update_nextclade)

    p_info = sub.add_parser("info", help="Show architecture and class names")
    p_info.set_defaults(func=_cmd_info)
    return parser


def main(argv: list[str] | None = None) -> None:
    parser = build_parser()
    args = parser.parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
