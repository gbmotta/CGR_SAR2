# CGR-SAR2

Identificacao de variantes de SARS-CoV-2 a partir do genoma completo.

Ha **dois modos**:

1. **Atual (recomendado)** — linhagem **Pango** e clado **Nextstrain** via [Nextclade](https://docs.nextstrain.org/projects/nextclade/), com o dataset SARS-CoV-2 mais recente. Painel de circulacao publica (CovSpectrum / GenBank) e estatuto WHO (VOI/VUM de julho de 2026).
2. **Historico (artigo)** — CNN + CGR do Experimento 1 (Alpha, Beta, Gamma, Delta, Iota, Epsilon; GISAID 2021). Nao serve para Omicron nem para as linhagens de 2025–2026.

## App Streamlit (conta share.streamlit.io)

Localmente:

```bash
python -m cgr_sar2 update-nextclade
pip install -r requirements-streamlit.txt
python -m cgr_sar2 streamlit
```

Na [Streamlit Community Cloud](https://share.streamlit.io):

1. Publique este repositorio no GitHub.
2. Em share.streamlit.io ? **New app** ? escolha o repo.
3. Main file: `streamlit_app.py`
4. Advanced ? Requirements file: `requirements-streamlit.txt`
5. O `setup.sh` instala o binario Nextclade e descarrega o dataset atual no build.

Para refrescar linhagens depois de um update da Pango: **Reboot** da app (volta a correr o `setup.sh`) ou, em local, `python -m cgr_sar2 update-nextclade`.

## Linha de comando (Pango atual)

```bash
python -m cgr_sar2 update-nextclade
python -m cgr_sar2 pango --fasta genoma.fasta
```

## Modo historico (artigo IEEE Access)

```bash
pip install -r requirements.txt
python -m cgr_sar2 train --from-original --epochs 40
python -m cgr_sar2 predict --fasta examples/delta_example.fasta
python -m cgr_sar2 app
```

O modelo em `models/variant_cnn.pt` foi treinado nas 42.638 imagens CGR do subset 1A (GISAID). Holdout 10%: 100% de acerto, em linha com o artigo. Esse classificador **nao** acompanha XFG, NB.1.8.1, JN.1, etc.

## Circulacao WHO (julho 2026)

| Pango | Nextstrain | WHO |
|---|---|---|
| JN.1 | 24A | VOI |
| XFG | 25C | VUM |
| NB.1.8.1 | 25B | VUM |
| PQ.16.1.1 | 25B | VUM |
| BA.3.2 | — | VUM |

Fonte: [WHO Tracking SARS-CoV-2 variants](https://www.who.int/activities/tracking-SARS-CoV-2-variants). A atribuicao de uma amostra concreta e feita pelo Nextclade (arvore de referencia atualizada), nao por esta tabela.

## Dados e scripts originais do artigo

Em `original/matlab`, `original/R`, `original/python` e `data/original/` (imagens `.mat` do Dropbox).

## Referencia

Camara, G. B. M., Coutinho, M. G. F., Barbosa, R. de M., & Fernandes, M. A. C. (2026). Integrating Chaos Game Representation with Convolutional Neural Networks for Accurate Classification of SARS-CoV-2 Variants. *IEEE Access*. DOI: [10.1109/ACCESS.2026.3723929](https://doi.org/10.1109/ACCESS.2026.3723929)
