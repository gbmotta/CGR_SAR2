"""Treino da CNN nas imagens CGR 64x64 (pipeline MATLAB original)."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np
import torch
from sklearn.metrics import classification_report, confusion_matrix
from sklearn.model_selection import train_test_split
from torch import nn
from torch.utils.data import DataLoader, Dataset

from .cgr import ambiguous_fraction, sequence_to_fcgr
from .fasta import parse_fasta
from .labels import (
    MAX_AMBIGUOUS_FRACTION,
    MIN_GENOME_LENGTH,
    VARIANT_LABELS,
    VARIANT_TO_INDEX,
    who_name_from_gisaid,
)
from .model import VariantCNN, count_parameters

DEFAULT_DATA_DIR = Path("data/fasta")
DEFAULT_CACHE_DIR = Path("data/fcgr")
DEFAULT_ORIGINAL_DIR = Path("data/original")
DEFAULT_OUT = Path("models/variant_cnn.pt")


class FcgrDataset(Dataset):
    def __init__(self, images: np.ndarray, labels: np.ndarray) -> None:
        self.images = torch.from_numpy(images).unsqueeze(1)
        self.labels = torch.from_numpy(labels.astype(np.int64))

    def __len__(self) -> int:
        return int(self.labels.shape[0])

    def __getitem__(self, index: int) -> tuple[torch.Tensor, torch.Tensor]:
        return self.images[index], self.labels[index]


def load_labeled_fastas(data_dir: Path) -> tuple[list[np.ndarray], list[int], list[str]]:
    images: list[np.ndarray] = []
    labels: list[int] = []
    ids: list[str] = []
    for variant in VARIANT_LABELS:
        folder = data_dir / variant.lower()
        if not folder.is_dir():
            continue
        for fasta_path in sorted(folder.glob("*.fa*")):
            for header, seq in parse_fasta(fasta_path):
                if len(seq) < MIN_GENOME_LENGTH:
                    continue
                if ambiguous_fraction(seq) > MAX_AMBIGUOUS_FRACTION:
                    continue
                images.append(sequence_to_fcgr(seq))
                labels.append(VARIANT_TO_INDEX[variant])
                ids.append(f"{variant}|{header.split()[0]}")
    if not images:
        raise FileNotFoundError(
            f"Nenhuma sequencia valida em {data_dir}/<variante>/*.fasta. "
            "Execute primeiro `python -m cgr_sar2 download` "
            "ou treine com --from-original."
        )
    return images, labels, ids


def load_or_build_cache(data_dir: Path, cache_dir: Path) -> tuple[np.ndarray, np.ndarray, list[str]]:
    cache_dir.mkdir(parents=True, exist_ok=True)
    x_path = cache_dir / "images.npy"
    y_path = cache_dir / "labels.npy"
    id_path = cache_dir / "ids.json"
    stamp = cache_dir / "source.txt"
    marker = f"{data_dir.resolve()}:{_fasta_count(data_dir)}:matlab-cgr"
    if x_path.is_file() and stamp.is_file() and stamp.read_text() == marker:
        return np.load(x_path), np.load(y_path), json.loads(id_path.read_text())
    images, labels, ids = load_labeled_fastas(data_dir)
    x = np.stack(images).astype(np.float32)
    y = np.asarray(labels, dtype=np.int64)
    np.save(x_path, x)
    np.save(y_path, y)
    id_path.write_text(json.dumps(ids))
    stamp.write_text(marker)
    return x, y, ids


def _fasta_count(data_dir: Path) -> int:
    return sum(1 for _ in data_dir.glob("*/*.fa*"))


def _mat_to_nchw(arr: np.ndarray) -> np.ndarray:
    """MATLAB (H, W[, 1], N) uint8 -> (N, H, W) float32 in [0, 1]."""
    if arr.ndim == 4:
        # (H, W, C, N)
        images = np.transpose(arr, (3, 0, 1, 2))[:, :, :, 0]
    elif arr.ndim == 3:
        # (H, W, N)
        images = np.transpose(arr, (2, 0, 1))
    else:
        raise ValueError(f"Unexpected CGR array shape: {arr.shape}")
    return images.astype(np.float32) / 255.0


def _labels_from_csv(path: Path) -> np.ndarray:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    names = [who_name_from_gisaid(row["Variant"]) for row in rows]
    return np.asarray([VARIANT_TO_INDEX[name] for name in names], dtype=np.int64)


def load_original_experiment1(
    original_dir: Path,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    import scipy.io as sio

    train_mat = original_dir / "val_images.mat"
    test_mat = original_dir / "teste_10pct_images.mat"
    train_csv = original_dir / "subset_1a_labels.csv"
    test_csv = original_dir / "subset_1b_labels.csv"
    missing = [p for p in (train_mat, test_mat, train_csv, test_csv) if not p.is_file()]
    if missing:
        raise FileNotFoundError(
            "Imagens originais em falta: " + ", ".join(str(p) for p in missing)
        )
    train_dict = sio.loadmat(train_mat)
    test_dict = sio.loadmat(test_mat)
    x_train = _mat_to_nchw(train_dict["seqsCGR_train"])
    x_test = _mat_to_nchw(test_dict["seqsCGR_teste"])
    y_train = _labels_from_csv(train_csv)
    y_test = _labels_from_csv(test_csv)
    if x_train.shape[0] != y_train.shape[0]:
        raise ValueError(
            f"Subset 1A: {x_train.shape[0]} imagens vs {y_train.shape[0]} rotulos"
        )
    if x_test.shape[0] != y_test.shape[0]:
        raise ValueError(
            f"Subset 1B: {x_test.shape[0]} imagens vs {y_test.shape[0]} rotulos"
        )
    return x_train, y_train, x_test, y_test


def train_model(
    data_dir: Path = DEFAULT_DATA_DIR,
    cache_dir: Path = DEFAULT_CACHE_DIR,
    out_path: Path = DEFAULT_OUT,
    epochs: int = 40,
    batch_size: int = 128,
    learning_rate: float = 1e-3,
    seed: int = 42,
    device: str | None = None,
    from_original: bool = False,
    original_dir: Path = DEFAULT_ORIGINAL_DIR,
    zscore: bool = False,
) -> Path:
    torch.manual_seed(seed)
    np.random.seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)
    device_t = torch.device(device or ("cuda" if torch.cuda.is_available() else "cpu"))
    pixel_mean = None
    pixel_std = None

    if from_original:
        x_pool, y_pool, _x_holdout, _y_holdout = load_original_experiment1(original_dir)
        print(
            f"Experimento 1 original: {x_pool.shape[0]} imagens alinhadas (subset 1A). "
            "O CSV do subset 1B nao corresponde a teste_10pct_images.mat; "
            "avaliacao num holdout estratificado de 10% de 1A (como um fold do k=10)."
        )
        x_train, x_val, y_train, y_val = train_test_split(
            x_pool, y_pool, test_size=0.10, random_state=seed, stratify=y_pool
        )
        x_test, y_test = x_val, y_val
        source = "GISAID Experiment 1 (Dropbox val_images.mat, holdout 10% de 1A)"
    else:
        x, y, _ids = load_or_build_cache(data_dir, cache_dir)
        x_train, x_test, y_train, y_test = train_test_split(
            x, y, test_size=0.10, random_state=seed, stratify=y
        )
        x_train, x_val, y_train, y_val = train_test_split(
            x_train, y_train, test_size=1 / 9, random_state=seed, stratify=y_train
        )
        source = f"FASTA em {data_dir}"

    counts = {VARIANT_LABELS[i]: int((y_train == i).sum()) for i in range(len(VARIANT_LABELS))}
    print("Amostras de treino por classe:", counts)
    print(f"Dispositivo: {device_t}")

    if zscore:
        pixel_mean = x_train.mean(axis=0, keepdims=True)
        pixel_std = np.maximum(x_train.std(axis=0, keepdims=True), 1e-6)
        x_train = (x_train - pixel_mean) / pixel_std
        x_val = (x_val - pixel_mean) / pixel_std
        x_test = (x_test - pixel_mean) / pixel_std

    train_loader = DataLoader(
        FcgrDataset(x_train, y_train),
        batch_size=batch_size,
        shuffle=True,
        pin_memory=device_t.type == "cuda",
    )
    val_loader = DataLoader(FcgrDataset(x_val, y_val), batch_size=batch_size)
    test_loader = DataLoader(FcgrDataset(x_test, y_test), batch_size=batch_size)

    model = VariantCNN().to(device_t)
    print(f"Parametros treinaveis: {count_parameters(model):,}")
    optimizer = torch.optim.Adam(model.parameters(), lr=learning_rate)
    criterion = nn.CrossEntropyLoss()

    best_acc = -1.0
    best_state = None
    out_path.parent.mkdir(parents=True, exist_ok=True)

    for epoch in range(1, epochs + 1):
        model.train()
        running = 0.0
        n_seen = 0
        for batch_x, batch_y in train_loader:
            batch_x = batch_x.to(device_t, non_blocking=True)
            batch_y = batch_y.to(device_t, non_blocking=True)
            optimizer.zero_grad()
            logits = model(batch_x)
            loss = criterion(logits, batch_y)
            loss.backward()
            optimizer.step()
            running += float(loss.item()) * batch_x.size(0)
            n_seen += batch_x.size(0)
        val_acc, val_loss = _evaluate(model, val_loader, criterion, device_t)
        print(
            f"Epoca {epoch:02d}/{epochs}  loss={running / n_seen:.4f}  "
            f"val_loss={val_loss:.4f}  val_acc={val_acc:.4f}"
        )
        if val_acc > best_acc:
            best_acc = val_acc
            best_state = {k: v.detach().cpu().clone() for k, v in model.state_dict().items()}

    assert best_state is not None
    model.load_state_dict(best_state)
    test_acc, _ = _evaluate(model, test_loader, criterion, device_t)
    y_true, y_pred = _collect(model, test_loader, device_t)
    print("\nConjunto de teste")
    print(classification_report(y_true, y_pred, target_names=list(VARIANT_LABELS), digits=4, zero_division=0))
    print("Matriz de confusao:\n", confusion_matrix(y_true, y_pred))
    payload = {
        "model": best_state,
        "labels": list(VARIANT_LABELS),
        "val_accuracy": best_acc,
        "test_accuracy": test_acc,
        "architecture": "MATLAB Experiment 1 VariantCNN (7x7/5x5/5x5)",
        "k": 6,
        "image_size": 64,
        "num_parameters": count_parameters(model),
        "source": source,
        "cgr": "MATLAB CGR.m + histcounts2 + rot90 + 8-bit",
    }
    if pixel_mean is not None and pixel_std is not None:
        payload["pixel_mean"] = pixel_mean.astype(np.float32)
        payload["pixel_std"] = pixel_std.astype(np.float32)
    torch.save(payload, out_path)
    print(f"Modelo guardado em {out_path} (val_acc={best_acc:.4f}, test_acc={test_acc:.4f})")
    return out_path


@torch.no_grad()
def _evaluate(
    model: VariantCNN,
    loader: DataLoader,
    criterion: nn.Module,
    device: torch.device,
) -> tuple[float, float]:
    model.eval()
    correct = 0
    total = 0
    loss_sum = 0.0
    for batch_x, batch_y in loader:
        batch_x = batch_x.to(device)
        batch_y = batch_y.to(device)
        logits = model(batch_x)
        loss_sum += float(criterion(logits, batch_y).item()) * batch_x.size(0)
        pred = logits.argmax(dim=1)
        correct += int((pred == batch_y).sum().item())
        total += batch_x.size(0)
    return correct / max(total, 1), loss_sum / max(total, 1)


@torch.no_grad()
def _collect(model: VariantCNN, loader: DataLoader, device: torch.device) -> tuple[list[int], list[int]]:
    model.eval()
    y_true: list[int] = []
    y_pred: list[int] = []
    for batch_x, batch_y in loader:
        logits = model(batch_x.to(device))
        y_true.extend(batch_y.tolist())
        y_pred.extend(logits.argmax(dim=1).cpu().tolist())
    return y_true, y_pred


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Treina a CNN CGR de variantes de SARS-CoV-2")
    parser.add_argument("--data-dir", type=Path, default=DEFAULT_DATA_DIR)
    parser.add_argument("--cache-dir", type=Path, default=DEFAULT_CACHE_DIR)
    parser.add_argument("--original-dir", type=Path, default=DEFAULT_ORIGINAL_DIR)
    parser.add_argument(
        "--from-original",
        action="store_true",
        help="Treinar nas imagens MATLAB do Experimento 1 (data/original/*.mat)",
    )
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--epochs", type=int, default=40)
    parser.add_argument("--batch-size", type=int, default=128)
    parser.add_argument("--lr", type=float, default=1e-3)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument(
        "--zscore",
        action="store_true",
        help="Normalizacao z-score por pixel (nao usada no MATLAB original)",
    )
    return parser


def main(argv: list[str] | None = None) -> None:
    args = build_parser().parse_args(argv)
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


if __name__ == "__main__":
    main()
