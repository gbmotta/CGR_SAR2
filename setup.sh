#!/usr/bin/env bash
# Streamlit Community Cloud corre este script no build.
set -euo pipefail
mkdir -p bin data/nextclade
ARCH="$(uname -m)"
case "$ARCH" in
  x86_64|amd64) ASSET="nextclade-x86_64-unknown-linux-gnu" ;;
  aarch64|arm64) ASSET="nextclade-aarch64-unknown-linux-gnu" ;;
  *) echo "Arquitectura nao suportada: $ARCH" >&2; exit 1 ;;
esac
curl -fsSL -o bin/nextclade \
  "https://github.com/nextstrain/nextclade/releases/latest/download/${ASSET}"
chmod +x bin/nextclade
./bin/nextclade dataset get \
  --name nextstrain/sars-cov-2/wuhan-hu-1/orfs \
  --output-dir data/nextclade/sars-cov-2
./bin/nextclade --version
