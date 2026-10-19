# Comorbidity Network NLP

> Code repository for:
> **"Enhancing comorbidity network inference with risk-enriched health trajectories embedding"**
> *(Fontana N. & Mapelli A., Di Angelantonio E., Ieva F.)*

This repository implements a pipeline for estimating **cardiometabolic comorbidity networks** from longitudinal patient event sequences, explicitly adjusted for confounding by shared clinical risk factors. Diseases and risk factors are embedded as dense vectors by fine-tuning a pretrained biomedical language model (BioWordVec) on patient trajectories. Pairwise cosine similarities between embeddings are used to detect disease-specific confounders, which are then incorporated as prior knowledge into a sparse Gaussian Graphical Model (GGM). The resulting network is analysed topologically, and patient-level disease progression dynamics are characterised through community transition matrices, clustering, and survival analysis.

Because the original patient data are confidential, a **synthetic sequence simulator** is provided that reproduces the structure and epidemiology of the real cohort.

---

## Repository structure

```
Com_network_NLP/
├── ComNetworkNLP.Rproj             ← Open this in RStudio for all R steps
├── environment.yml                 ← Conda environment (Python steps)
├── requirements.txt                ← Pip alternative
│
├── simulate_sequences.py           ← Generate synthetic patient sequences
│
├── BioWordVec/
│   ├── BioWordVec.py               ← Train base BioWordVec model (run once)
│   └── data/
│       ├── MeSH_dic.pkl.gz         ← MeSH ID → term dictionary
│       ├── MeSH_graph.edgelist     ← MeSH taxonomy graph
│       └── pubmed_sample           ← PubMed abstracts corpus
│
├── 1.Sequence_preprocessing.ipynb  ← Step 1: sequence construction (Python)
├── 2.Word_embedding.ipynb          ← Step 2: BioWordVec fine-tuning (Python)
├── 3.Prior_clinical_knowledge.ipynb← Step 3: confounder detection (Python)
├── 4.Network_estimation.R          ← Step 4: GGM network estimation (R)
├── 5.Topological_analysis.Rmd      ← Step 5: topology & communities (R)
├── Downstream_analysis.Rmd         ← Patient clustering & survival (R)
│
├── functions/
│   └── ggReg_main_functions.R      ← GGM neighbourhood selection functions
│
├── data_simulated/                 ← Created by simulate_sequences.py
├── results/                        ← Shared intermediate outputs across steps
└── plots/                          ← All generated figures
```

---

## Pipeline overview

| Step | File | Language | What it does |
|---|---|---|---|
| 0 | `simulate_sequences.py` | Python | Generate synthetic patient sequences |
| — | `BioWordVec/BioWordVec.py` | Python | Pretrain base BioWordVec (run once) |
| 1 | `1.Sequence_preprocessing.ipynb` | Python | Build chronologically ordered disease and risk-factor sequences; filter cohort |
| 2 | `2.Word_embedding.ipynb` | Python | Fine-tune BioWordVec on patient sequences; extract disease and risk-factor embeddings |
| 3 | `3.Prior_clinical_knowledge.ipynb` | Python | Compute differential similarity matrix; detect and validate confounders; build prior knowledge matrix |
| 4 | `4.Network_estimation.R` | R | Estimate sparse GGM network using Lasso neighbourhood selection with confounder priors; plot binary adjacency matrix; compute and plot the percentage of patients exhibiting each connected disease transition |
| 5 | `5.Topological_analysis.Rmd` | R | Compute centrality measures; detect disease communities via modularity optimisation |
| — | `Downstream_analysis.Rmd` | R | Build community transition matrices; cluster patients via k-means; Kaplan–Meier survival analysis |

Each Python step reads from and writes to the shared `results/` folder. Figures are saved to `plots/`. The R steps load their inputs from `results/` as well.

---

## Quick start (simulated data)

### 1. Set up the Python environment

```bash
conda env create -f environment.yml
conda activate comorb-network
```

### 2. Generate synthetic patient sequences (~30 seconds)

```bash
python simulate_sequences.py --n_patients 5000 --seed 42
```

This creates three files in `data_simulated/`:
- `sequence_main_simulated.csv` — disease events
- `sequence_parallel_simulated.csv` — risk-factor events
- `survival_outcome_simulated.csv` — death/censoring outcome per patient

### 3. Pretrain the BioWordVec base model (~15–30 min, run once)

```bash
cd BioWordVec
python BioWordVec.py
cd ..
```

This saves `pubmed_mesh_test` and `pubmed_mesh_test.bin` inside `BioWordVec/`. These files are required by step 2 and only need to be generated once.

### 4. Run the Python pipeline steps in order

Open and execute the three notebooks sequentially:

```
1.Sequence_preprocessing.ipynb
2.Word_embedding.ipynb
3.Prior_clinical_knowledge.ipynb
```

### 5. Run the R pipeline steps

> **Always open `ComNetworkNLP.Rproj` in RStudio first.** This sets the working directory to the repo root so all relative paths in the R files resolve correctly. Do not open the R or Rmd files directly by double-clicking.

Run in order:

| File | How to run in RStudio |
|---|---|
| `4.Network_estimation.R` | Source the script (`Ctrl+Shift+S`) |
| `5.Topological_analysis.Rmd` | Knit to HTML (`Ctrl+Shift+K`) |
| `Downstream_analysis.Rmd` | Knit to HTML (`Ctrl+Shift+K`) |

---

## Using real data

If you have access to real patient sequences, skip `simulate_sequences.py` and configure the data paths at the top of `1.Sequence_preprocessing.ipynb`:

```python
DATA_MODE              = "real"
REAL_MAIN_SEQ_PATH     = "/path/to/your/sequence_main.csv"
REAL_PARALLEL_SEQ_PATH = "/path/to/your/sequence_parallel.csv"
```

### Required input format

Both files must have exactly these three columns:

| Column | Type | Description |
|---|---|---|
| `ID` | int or str | Patient identifier |
| `age_in_days_since_born` | int | Age at event occurrence (days since birth) |
| `event_name` | str | MeSH disease or risk-factor name |

**Main sequence** (`sequence_main.csv`): one row per disease event, covering the 24 cardiometabolic MeSH disease names (e.g. `Hypertension`, `Coronary Disease`, `Heart Failure`). Events must be temporally ordered within each patient.

**Parallel sequence** (`sequence_parallel.csv`): one row per risk-factor observation, covering demographics, biomarkers, comorbidities, and socioeconomic variables (e.g. `Diabetes Control, Good`, `Body Mass Index, High`, `Hypertension Control, Bad`).

The simulated files in `data_simulated/` are the best reference for the expected vocabulary and format.

---

## User-configurable settings

Each script has a configuration section that should be reviewed before running. Key decision points are highlighted below.

### Step 1 — `1.Sequence_preprocessing.ipynb` (top cell)

```python
DATA_MODE        = "simulated"   # "simulated" | "real"
OUTPUT_RES_DIR   = "results"     # folder for all intermediate outputs
OUTPUT_PLOTS_DIR = "plots"       # folder for all figures
```

### Step 2 — `2.Word_embedding.ipynb`

```python
FINETUNE_EPOCHS = 4
# Number of fine-tuning epochs applied to BioWordVec.
# The paper uses 4, selected so that mean cosine similarity
# between original and fine-tuned embeddings stays above 0.70.
# More epochs → stronger sequence adaptation, less biomedical prior retained.
```

### Step 3 — `3.Prior_clinical_knowledge.ipynb`

```python
RF_CORRELATION_THRESHOLD = <median of RF–RF similarity submatrix>
# A risk factor is considered associated with a disease if their cosine
# similarity exceeds this threshold. Default is the median pairwise
# similarity among all risk factors (mirrors the paper criterion).

TOP_N_PAIRS = <95th percentile cutoff>
# Disease pairs with |Δ similarity| above this cutoff are treated as
# "highly affected" and analysed for shared confounders.
# Default is all pairs above the 95th percentile of absolute changes.
```

### Step 4 — `4.Network_estimation.R` (top configuration block)

```r
OUTPUT_RES_DIR  = "results"    # folder for intermediate outputs (RData, CSV)
OUTPUT_PLOT_DIR = "plots"      # folder for figures
color_net_plots = "#FFA420"    # fill colour used in the adjacency matrix heatmap
```

This script produces three outputs in addition to the GGM model itself:
- `plots/disease_adj_matrix_GGM_heatmap.svg/.png` — binary lower-triangle heatmap showing which disease pairs have an estimated edge
- `plots/disease_adj_matrix_GGM_heatmap_transition_percentage.svg/.png` — heatmap showing, for each connected disease pair, the percentage of patients who exhibit at least one consecutive transition between the two diseases in their event sequence
- `results/GGM_transition_analysis.RData` — R objects `adj_GGM_df`, `transition_counts`, and `transition_percentage` for downstream use

### Step 5 — `5.Topological_analysis.Rmd` (YAML params block at the top of the file)

```yaml
params:
  OUTPUT_RES_DIR: "results"
  OUTPUT_PLOT_DIR: "plots"
  community_detection_algorithm: "cfg_weighted"
  # Community detection algorithm. Options:
  #   cfg           — fast greedy modularity optimisation (unweighted)
  #   cfg_weighted  — fast greedy modularity optimisation (weighted)  ← paper default
  #   leid          — Leiden multi-level optimisation (unweighted)
  #   leid_weighted — Leiden multi-level optimisation (weighted)
  #   btw_cluster   — edge betweenness clustering
  color_net_plots: "#FFA420"   # node colour in network plots
```

### Downstream analysis — `Downstream_analysis.Rmd` (YAML params block + user input cell)

```yaml
params:
  OUTPUT_RES_DIR:    "results"
  OUTPUT_PLOT_DIR:   "plots"
  OUTCOME_FILE_PATH: "data_simulated/survival_outcome_simulated.csv"
  # Survival outcome file — must have columns: ID, event_time, event_status
  # event_status: 1 = death, 0 = censored
  # event_time:   age in days at death or censoring
  MAPPING_COM_PATH: "sim_mapping_communities.csv"
  # Optional: a CSV with columns community_number and community_name
  # providing human-readable labels for the detected communities.
  # Create this file manually after inspecting the community assignments
  # produced by step 5. Leave as "" to use auto-labels (COM_1, COM_2, …).
```

After the elbow plot is rendered, you must set the number of patient clusters:

> **User input required:** inspect the elbow plot and set `k_chosen` in the following cell.

```r
k_chosen = 4   # replace with your chosen number of clusters
```

---

## Dependencies

### Python

```bash
conda env create -f environment.yml   # recommended
# or
pip install -r requirements.txt
```

| Package | Version | Purpose |
|---|---|---|
| `gensim` | ≥ 4.3 | FastText / BioWordVec embeddings |
| `node2vec` | ≥ 0.4 | MeSH graph random walks (BioWordVec pretraining) |
| `networkx` | ≥ 3.0 | Graph utilities |
| `scikit-learn` | ≥ 1.3 | Cosine similarity matrices |
| `scipy` | ≥ 1.10 | Statistical tests (Wilcoxon, Mann–Whitney, KS) |
| `matplotlib`, `seaborn` | — | Visualisation |
| `leidenalg`, `igraph` | ≥ 0.10 | Leiden community detection (optional) |
| `python-louvain` | ≥ 0.16 | Louvain community detection (optional) |

### R

Install the required packages from CRAN before running steps 4–5 and the downstream analysis:

```r
install.packages(c(
  # Network analysis & visualisation
  "igraph", "tidygraph", "ggraph", "backbone", "poweRlaw",
  # Data wrangling
  "dplyr", "tidyr", "tibble", "readr",
  # Plotting
  "ggplot2", "RColorBrewer", "patchwork", "scales", "gridExtra", "plotly",
  # Tables
  "knitr", "kableExtra",
  # Survival & clustering
  "survival", "survminer", "cluster"
))
```

---

## Simulated data disclaimer

> The synthetic sequences are for **demonstration only.**
> Simulated patients do not represent real individuals, and network structures produced from simulated data will differ from those reported in the original publication.

The simulator reproduces the real cohort's:
- 24 cardiometabolic MeSH disease types with empirical incidence rates (e.g. Hypertension 96%, Dyslipidemias 75%)
- ~81–97 parallel risk-factor types spanning demographics, biomarkers, comorbidities, and socioeconomic factors
- Temporal ordering constraints (e.g. CKD Stage 3 → Stage 4 → Stage 5 → ESRD)
- Right-censored survival outcome driven by disease severity score

---

## Citation

If you use this code, please cite:

```bibtex
@article{fontana2026enhancing,
  title={Enhancing comorbidity network inference with risk-enriched health trajectories embedding},
  author={Fontana, Nicole and Mapelli, Alessia and Di Angelantonio, Emanuele and Ieva, Francesca},
  journal={arXiv preprint arXiv:2607.04702},
  year={2026}
}
```

---

## License

MIT License — see `LICENSE` for details.
