# Differential Gene Expression Analysis and Visualization Pipeline for Cancer Data

## 📋 Description

This project is a streamlined **bioinformatics pipeline** designed to perform **differential gene expression analysis** and visual exploration of gene relationships in cancer datasets.

The core of the project is built around the `pydeseq2` library, which enables statistical testing for **differentially expressed genes (DEGs)** using DESeq2 methodology — a standard in RNA-Seq data analysis. 

The project takes as input:
- A **gene expression count matrix**.
- Corresponding **sample metadata**.

It performs automated DESeq2-based analysis and generates ranked results based on **adjusted p-values**.

---

## 🎨 Visualizations Included

In addition to statistical analysis, the project features an `OmicsVisualizer` class offering a suite of publication-ready plots for exploratory data analysis:

- 📈 **Volcano Plot**  
  Highlights **upregulated** and **downregulated** genes based on statistical significance and fold-change thresholds.

- 🔥 **Heatmap**  
  Visualizes expression patterns for the **top 25 differentially expressed genes** across patient samples.

- 🔗 **Correlation Matrix Heatmap**  
  Displays **gene-gene correlation** to reveal potential co-expression patterns.

- 🧬 **Circular Plot**  
  Offers a **network-based view** of gene relationships filtered by user-defined correlation cutoffs.
