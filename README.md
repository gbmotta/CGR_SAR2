# CGR-SAR2

Identificacao de variantes de SARS-CoV-2 a partir do genoma completo.

Ha **dois modos**:

1. **Atual (recomendado)** — linhagem **Pango** e clado **Nextstrain** via [Nextclade](https://docs.nextstrain.org/projects/nextclade/). Opcionalmente **Pangolin/UShER** (Docker). Painel de circulacao (CovSpectrum) e mutacoes-assinatura WHO (XFG, NB.1.8.1, PQ.16.1.1, JN.1).
2. **Historico (artigo)** — CNN + CGR do Experimento 1 (Alpha–Epsilon, GISAID 2021). Nao serve para Omicron nem para 2025–2026.

## App Streamlit

```bash
python -m cgr_sar2 update-nextclade
pip install -r requirements-streamlit.txt
python -m cgr_sar2 streamlit
```

Aceita **varios FASTA**, mostra assinaturas da espicula e exporta **CSV/Excel**.

Na [Streamlit Community Cloud](https://share.streamlit.io): repo `gbmotta/CGR_SAR2`, main file `streamlit_app.py`, requirements `requirements-streamlit.txt`. O `setup.sh` instala o Nextclade no build.

## Linha de comando (lote)

```bash
python -m cgr_sar2 update-nextclade
python -m cgr_sar2 pango --fasta amostra1.fasta amostra2.fasta --csv saida.csv --xlsx saida.xlsx
python -m cgr_sar2 pango --input-dir pasta_com_fastas --engine both
```

`--engine both` corre Nextclade e Pangolin (se o binario ou `staphb/pangolin` estiver disponivel) e indica se as linhagens coincidem.

## API

```bash
pip install -r requirements-api.txt
python -m cgr_sar2 api --port 8000
```

```bash
curl -s http://127.0.0.1:8000/health
curl -s -X POST http://127.0.0.1:8000/classify \
  -H 'Content-Type: application/json' \
  -d '{"fasta":">q\nATTAAAGGTTTATACCTTCCCAGGTAACAA...","engine":"nextclade"}'
curl -s -X POST -F fasta=@genoma.fasta http://127.0.0.1:8000/classify/file
```

## Docker

```bash
docker compose up --build api
# UI: docker compose up --build ui
# Pangolin pontual:
docker compose --profile lab run --rm pangolin /data/query.fasta --outfile /data/lineage_report.csv
```

A imagem inclui Nextclade. Pangolin e um servico opcional (`staphb/pangolin`) porque o UShER e pesado demais para a Cloud.

## Dataset Nextclade (semanal)

O workflow `.github/workflows/update-nextclade.yml` corre **todas as segundas** (ou manualmente) e faz commit do dataset se a Nextstrain publicou uma versao nova. Localmente: `python -m cgr_sar2 update-nextclade`.

## Exemplos atuais

```bash
python -m cgr_sar2 download --current-examples
```

Gera `examples/xfg_example.fasta`, `nb181_example.fasta` e `jn1_example.fasta` a partir do NCBI. Os FASTA Alpha–Epsilon ficam em `examples/historical/`.

## Modo historico (artigo IEEE Access)

```bash
pip install -r requirements.txt
python -m cgr_sar2 train --from-original --epochs 40
python -m cgr_sar2 predict --fasta examples/historical/delta_example.fasta
```

## Circulacao WHO (julho 2026)

| Pango | Nextstrain | WHO |
|---|---|---|
| JN.1 | 24A | VOI |
| XFG | 25C | VUM |
| NB.1.8.1 | 25B | VUM |
| PQ.16.1.1 | 25B | VUM |
| BA.3.2 | — | VUM |

## Referencia

Camara, G. B. M., Coutinho, M. G. F., Barbosa, R. de M., & Fernandes, M. A. C. (2026). Integrating Chaos Game Representation with Convolutional Neural Networks for Accurate Classification of SARS-CoV-2 Variants. *IEEE Access*. DOI: [10.1109/ACCESS.2026.3723929](https://doi.org/10.1109/ACCESS.2026.3723929)
