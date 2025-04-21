# 🧬 Gene Regulatory Network (GRN) Construction & Analysis

This project focuses on constructing a **Gene Regulatory Network** from gene expression data using known **Human Transcription Factors (TFs)** from [HumanTFs](https://humantfs.ccbr.utoronto.ca/). The GRN is modeled and analyzed using graph theory to uncover potential regulatory interactions between genes.

---

## 🌐 Dataset

- **Gene Expression Data**: A `.csv` file containing expression levels of genes.
- **Transcription Factors List**: Downloaded from [HumanTFs Database](https://humantfs.ccbr.utoronto.ca/download.php), which provides a curated list of human TFs and co-TFs.

---

## 🧪 Methodology

1. **Data Preparation**
   - Load gene expression data.
   - Normalize and clean data if necessary.
   - Map genes with known transcription factors using HumanTFs.

2. **Interaction Inference**
   - Compute pairwise Pearson correlation between transcription factors and target genes.
   - Apply thresholding to determine significant regulatory links.
   - Optionally use mutual information or partial correlation for robustness.

3. **Graph Construction**
   - Construct a **directed graph** using NetworkX where:
     - Nodes = Genes (TFs and targets)
     - Edges = Regulatory interactions (from TFs to target genes)

4. **Analysis & Visualization**
   - Visualize the network using NetworkX and Matplotlib.
   - Perform centrality analysis to identify key regulators.
   - Export graph metrics (e.g., degree, betweenness) for biological interpretation.

---

## 🛠️ Installation

Instructions to install and run the project.

Netowork.py is the main file for constructing the gene regulatory network.
For using it use the following commands
```bash
from Network import *

