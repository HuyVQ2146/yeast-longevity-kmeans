# Determining Potential Yeast Longevity Genes via Optimized K-Means Clustering

**Team Project** - **Team Members:** Vu Quang Huy, Pham Van Bach, Nguyen Nam Binh, Nguyen Tuan Minh

---

## Overview

This project applies **K-Means clustering** to yeast microarray gene expression data in order to identify **candidate longevity genes** — genes whose expression profiles are similar to already-known longevity genes and are therefore strong targets for future wet-lab validation.

The central idea: if a candidate gene clusters tightly with a high proportion of known longevity genes, it likely shares the same biological function and may also influence lifespan.

---

## Dataset

| File | Description |
|---|---|
| `Longotor1delta_test.xlsx` | Multi-sheet Excel workbook containing raw microarray data for ~5,667 *S. cerevisiae* genes across three conditions |
| `cleaned_data.csv` | Cleaned subset of 170 genes (87 candidate + 83 known longevity) matched from the raw data |
| `(READ ME!) interested_genes.txt` | Lists of candidate and known longevity gene names used to filter from the full dataset |

### Raw Data (`Longotor1delta_test.xlsx`)

The workbook contains four sheets:

| Sheet | Contents |
|---|---|
| **T. S2 microarray results** | Full dataset — 5,667 genes with expression log-ratios |
| **Count** | Helper sheet with gene-name counts |
| **Check** | Validation sheet cross-referencing candidate/known genes against the full dataset |
| **CLeaned Data** | The 170 matched genes with formulas referencing the raw sheet |

### Gene Expression Features

Each gene has three log-ratio expression values:

| Column | Description |
|---|---|
| `sch9/wt` | Expression ratio of **sch9Δ** mutant vs. wild-type |
| `ras2/wt` | Expression ratio of **ras2Δ** mutant vs. wild-type |
| `tor1/wt` | Expression ratio of **tor1Δ** mutant vs. wild-type |

These three pathways (Sch9, Ras2, Tor1) are well-characterized nutrient-sensing / longevity-regulating pathways in yeast.

---

## Methodology

### 1. Feature Engineering

- **Min-Max Normalization:** Each of the three expression features is rescaled from `[-10, 10]` → `[0, 1]`.
- **Euclidean Distance:** A single composite feature is computed as the Euclidean distance from the origin in the 3D normalized space:

  ```
  distance = √(norm_sch9² + norm_ras2² + norm_tor1²)
  ```

### 2. K-Means Clustering

- Standard Lloyd's algorithm with 100 max iterations.
- Initial centroids sampled randomly (seed = 42 for reproducibility).
- Applied for **K = 2 to 10**, using 1D distance as the clustering feature.

### 3. Optimal K Selection

- **Silhouette Score** is used as the evaluation metric.
- **K = 8** achieves the highest Silhouette Score of **0.5829**.

### 4. Cluster Analysis

Each cluster is profiled by the ratio of known longevity genes to total genes. The cluster with the highest density of known longevity genes is identified as the most promising group, and its candidate genes are flagged for experimental follow-up.

---

## Results

### Silhouette Scores

| K | Silhouette Score | Cluster Sizes |
|---|---|---|
| 2 | 0.5785 | [101, 69] |
| 3 | 0.5283 | [74, 56, 40] |
| 4 | 0.5422 | [68, 54, 13, 35] |
| 5 | 0.5643 | [55, 46, 12, 35, 22] |
| 6 | 0.5110 | [49, 41, 12, 32, 17, 19] |
| 7 | 0.5659 | [30, 41, 10, 33, 17, 10, 29] |
| **8** | **0.5829** | **[30, 41, 10, 32, 14, 7, 15, 21]** |
| 9 | 0.5503 | [33, 24, 19, 32, 14, 7, 15, 21, 5] |
| 10 | 0.5723 | [40, 32, 9, 13, 14, 7, 14, 10, 5, 26] |

### Cluster Distribution (K = 8)

| Cluster | Total Genes | Known Longevity | Candidate | Known Ratio |
|---|---|---|---|---|
| 1 | 30 | 10 | 20 | 33.3% |
| **2** | **41** | **27** | **14** | **65.9%** |
| 3 | 10 | 5 | 5 | 50.0% |
| 4 | 32 | 13 | 19 | 40.6% |
| 5 | 14 | 8 | 6 | 57.1% |
| 6 | 7 | 4 | 3 | 57.1% |
| 7 | 15 | 7 | 8 | 46.7% |
| 8 | 21 | 9 | 12 | 42.9% |
| **Total** | **170** | **83** | **87** | **48.8%** |

### Top Candidate Longevity Genes (Cluster 2)

Cluster 2 has the highest concentration of known longevity genes at **65.9%** (vs. the dataset average of 48.8%), making its 14 candidate genes the **strongest candidates** for new longevity genes:

| # | Gene ID / Name |
|---|---|
| 1 | YPL254W / HFI1 |
| 2 | YDR264C / AKR1 |
| 3 | YDR456W / NHX1 |
| 4 | YMR161W / HLJ1 |
| 5 | YBL007C / SLA1 |
| 6 | YMR297W / PRC1 |
| 7 | YOR209C / NPT1 |
| 8 | YDR069C / DOA4 |
| 9 | YOR023C / AHC1 |
| 10 | YAL020C / ATS1 |
| 11 | YBL047C / EDE1 |
| 12 | YKR052C / MRS4 |
| 13 | YKL205W / LOS1 |
| 14 | YGR200C / ELP2 |

---

## How to Run

### Prerequisites

- Python 3.x
- pandas

### Steps

```bash
# Install dependencies
pip install pandas

# Run the analysis
python Code_final.py
```

The script reads `cleaned_data.csv` from the same directory and prints Silhouette scores, cluster distributions, and gene lists to the console.

---

## Project Structure

```
.
├── Code_final.py                    # Main analysis script
├── cleaned_data.csv                 # Cleaned gene expression dataset (170 genes)
├── Longotor1delta_test.xlsx         # Raw multi-sheet microarray data
├── (READ ME!) interested_genes.txt  # Lists of candidate & known longevity genes
├── report.docx                      # Final written report
└── README.md                        # This file
```

---

## References

- **S. cerevisiae** longevity pathways: Sch9 (S6K homolog), Ras2 (Ras-cAMP-PKA), Tor1 (mTOR homolog) — major nutrient-sensing regulators of aging in yeast.
- Silhouette Score: Rousseeuw, P. J. (1987). *Silhouettes: A graphical aid to the interpretation and validation of cluster analysis.* Journal of Computational and Applied Mathematics, 20, 53–65.
- K-Means: MacQueen, J. (1967). *Some methods for classification and analysis of multivariate observations.*
