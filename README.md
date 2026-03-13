# Comorbidity Network Analysis Pipeline

> Code repository for:
> **"Enhancing comorbidity network inference with risk-enriched health trajectories embedding"**
> *(Fontana N. & Mapelli A., Di Angelantonio E., Ieva F.)*

This repository provides a fully reproducible Python pipeline for constructing and analysing **cardiometabolic comorbidity networks** from patient event sequences. The approach uses:

1. **BioWordVec** pre-training (Node2Vec on MeSH graph + FastText on biomedical text)
2. **FastText fine-tuning** on patient disease/risk-factor sequences
3. **Cosine similarity matrices** between disease embeddings
4. **Graphical Lasso (GGM)** for sparse network estimation
5. **Community detection** and centrality analysis

Because the original patient data are confidential, a **synthetic sequence simulator** is provided that faithfully mimics the structure and epidemiology of the real cohort.

---

## Repository Structure

```
comorbidity-network-analysis/
├── README.md
├── requirements.txt
├── environment.yml
│
├── simulate_sequences.py           ← Generate synthetic patient sequences
├── BioWordVec.py                   ← Train BioWordVec base model (run once)
├── Comorbidity_Network_Analysis.ipynb  ← Main analysis notebook (all sections)
│
├── data/                           ← Required input files for BioWordVec
│   ├── MeSH_dic.pkl.gz             ← MeSH ID → term dictionary
│   ├── MeSH_graph.edgelist         ← MeSH taxonomy graph
│   ├── pubmed_sample               ← Biomedical text corpus (PubMed abstracts)
│   └── MeSHFeatureGeneratedByDeepWalk.csv
│
├── data_simulated/                 ← Output of simulate_sequences.py
│   ├── sequence_main_simulated.csv
│   ├── sequence_parallel_simulated.csv
│   └── simulation_summary.txt
│
├── results_GGM/                    ← Output of Section 5 (GGM)
│   └── GGM_adj_main_par_full.csv
│
├── plots/                          ← All generated figures (SVG + PDF)
└── Baseline/                       ← Community detection results CSV
```

---

## Quick Start (Simulated Data)

```bash
# 1. Create and activate environment
conda env create -f environment.yml
conda activate comorb-network

# 2. Generate synthetic sequences (≈ 30 seconds)
python simulate_sequences.py --n_patients 5000 --seed 42

# 3. Train the BioWordVec base model (≈ 15–30 min, run ONCE)
python BioWordVec.py

# 4. Open and run the main analysis notebook
jupyter notebook Comorbidity_Network_Analysis.ipynb
```

After step 3 you will find `pubmed_mesh_test` and `pubmed_mesh_test.bin` in the current directory. These are required by Section 2 of the notebook.

---

## Using Real Data

If you have access to real patient sequences in the same format, skip step 2 above and set in the notebook configuration cell:

```python
DATA_MODE              = "real"
REAL_MAIN_SEQ_PATH     = "/path/to/your/sequence_main.csv"
REAL_PARALLEL_SEQ_PATH = "/path/to/your/sequence_parallel.csv"
```

### Required column format (both files)

| Column | Type | Description |
|---|---|---|
| `ID` | int/str | Patient identifier |
| `age_in_days_since_born` | int | Age at event (days) |
| `event_name` | str | **MeSH** disease or risk-factor name |

The 24 expected main-disease MeSH names and the mapping from original event types are documented in `Comorb_network_seq_preprocessing.ipynb` (chunks 21 and 23).

---

## Notebook Sections

| Section | Description | User Decision Points |
|---|---|---|
| 0 | Configuration | `DATA_MODE` |
| 1 | Data import, filtering, MeSH mapping, cohort summary, visualisations | — |
| 2 | BioWordVec fine-tuning, phrase-embedding matrix | `FINETUNE_EPOCHS`, `FINETUNE_WINDOW` |
| 3 | 24×24 and 121×121 cosine-similarity matrices, heatmaps, network graphs | — |
| 4 | Difference matrices, risk-factor analysis, sensitivity analysis, statistical tests | `RF_CORRELATION_THRESHOLD`, `TOP_N_PAIRS`, `THRESHOLD_VALUES` |
| 5 | Graphical Lasso (GGM) — fully Python, no R required | — |
| 6 | GGM network visualisation (heatmaps + graphs) | — |
| 7 | Topological analysis, community detection, centrality | — |

---

## GGM Estimation (Section 5)

Network sparsification uses **scikit-learn `GraphicalLassoCV`** (5-fold cross-validation to select the regularisation parameter `alpha`). No R installation is required.

The regularisation parameter is selected automatically by cross-validation and printed in the cell output. Results are saved to `results_GGM/GGM_adj_main_par_full.csv`.

---

## Simulated Data Disclaimer

> **The synthetic sequences are for demonstration only.**
> Simulated patients do not represent any real individuals. Network structures produced from simulated data will differ from the results reported in the original publication.

The simulation reproduces the real data's:
- 24 cardiometabolic disease types (MeSH names)
- ~81–97 parallel risk-factor types (demographics, biomarkers, comorbidities, socioeconomic)
- Real-data incidence rates (Hypertension 96%, Dyslipidemias 75%, …)
- Temporal ordering constraints (e.g. CKD Stage 3 → 4 → 5 → ESRD)

---

## Dependencies

Install via conda (recommended):

```bash
conda env create -f environment.yml
```

Or via pip:

```bash
pip install -r requirements.txt
```

### Key Python packages

| Package | Version | Purpose |
|---|---|---|
| `gensim` | ≥ 4.3 | FastText word embeddings |
| `node2vec` | ≥ 0.4 | MeSH graph random walks (used by BioWordVec.py) |
| `networkx` | ≥ 3.0 | Graph construction and analysis |
| `scikit-learn` | ≥ 1.3 | GraphicalLassoCV (GGM) |
| `scipy` | ≥ 1.10 | Statistical tests |
| `matplotlib`, `seaborn` | — | Visualisation |
| `leidenalg`, `igraph` | — | Leiden community detection (optional) |
| `python-louvain` | — | Alternative community detection |

> **Note:** `leidenalg` and `igraph` are optional. If not installed, only Fast Greedy and Girvan-Newman community detection are used.

---

## Citation

If you use this code, please cite:

```bibtex
@article{}
```

---

## License

MIT License — see `LICENSE` file for details.
