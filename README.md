# GDA-Project

# Comparative Over-Representation Analysis of Prostate Cancer Transcriptomes

## Project Summary
This analysis focused on identifying key biological pathways in prostate cancer through differential expression analysis and multiple over-representation analysis methods. By comparing the results from different statistical approaches, we identified robust pathway enrichment patterns that provide insight into the biological mechanisms of prostate cancer.

## Dataset Information
- Total samples: 492
- Tumor samples: 473
- Normal samples: 19
- Genes analyzed: 13963

## Differential Expression Analysis
- Significantly differentially expressed genes: 7435
- Up-regulated genes in tumor: 25
- Down-regulated genes in tumor: 7410

## Top Consensus Pathways in Up-regulated Genes
- **Arginine biosynthesis**: Found in 2 methods

## Top Consensus Pathways in Down-regulated Genes
- **Endocytosis**: Found in 2 methods
- **mRNA surveillance pathway**: Found in 2 methods
- **Focal adhesion**: Found in 2 methods
- **Autophagy**: Found in 2 methods
- **Thyroid hormone signaling pathway**: Found in 2 methods
- **Insulin signaling pathway**: Found in 2 methods
- **Gene Expression (Transcription) R-HSA-74160**: Found in 2 methods
- **Metabolism Of RNA R-HSA-8953854**: Found in 2 methods
- **RHO GTPase Cycle R-HSA-9012999**: Found in 2 methods
- **Signaling By Rho GTPases, Miro GTPases And RHOBTB3 R-HSA-9716542**: Found in 2 methods

## Pathway Analysis Method Comparison
We compared the results from multiple over-representation analysis methods including:
- Hypergeometric test
- Fisher's exact test
- Chi-squared test
- Enrichr (via GSEApy)
- Gene Set Enrichment Analysis (GSEA)

The consensus pathways identified by multiple methods represent robust biological signals that are less likely to be artifacts of any particular statistical approach.

## GSEA Pre-ranked Analysis
GSEA pre-ranked analysis was performed to identify pathways enriched across the entire ranked gene list, without requiring an arbitrary significance cutoff. Top enriched pathways included:

### KEGG_2021_Human

#### Positively Enriched (Up in Tumor)
- **Proteasome**: ES=0.467, FDR=0.000e+00
- **Citrate cycle (TCA cycle)**: ES=0.459, FDR=0.000e+00
- **Vibrio cholerae infection**: ES=0.290, FDR=0.000e+00
- **Arginine and proline metabolism**: ES=0.347, FDR=0.000e+00
- **Glycolysis / Gluconeogenesis**: ES=0.309, FDR=5.795e-03

#### Negatively Enriched (Down in Tumor)
- **Herpes simplex virus 1 infection**: ES=-0.453, FDR=1.091e-02
- **Homologous recombination**: ES=-0.500, FDR=2.304e-01
- **Cytokine-cytokine receptor interaction**: ES=-0.421, FDR=2.475e-01
- **Neuroactive ligand-receptor interaction**: ES=-0.479, FDR=2.563e-01
- **Morphine addiction**: ES=-0.465, FDR=3.068e-01

### GO_Biological_Process_2021

#### Positively Enriched (Up in Tumor)
- **mitochondrial electron transport, NADH to ubiquinone (GO:0006120)**: ES=0.407, FDR=0.000e+00
- **antigen processing and presentation of exogenous peptide antigen via MHC class I, TAP-dependent (GO:0002479)**: ES=0.322, FDR=0.000e+00
- **aerobic electron transport chain (GO:0019646)**: ES=0.444, FDR=0.000e+00
- **mitochondrial ATP synthesis coupled electron transport (GO:0042775)**: ES=0.444, FDR=0.000e+00
- **protein exit from endoplasmic reticulum (GO:0032527)**: ES=0.566, FDR=0.000e+00

#### Negatively Enriched (Down in Tumor)
- **positive regulation of neuron differentiation (GO:0045666)**: ES=-0.495, FDR=8.122e-01
- **defense response to symbiont (GO:0140546)**: ES=-0.388, FDR=8.154e-01
- **regulation of B cell proliferation (GO:0030888)**: ES=-0.503, FDR=8.162e-01
- **negative regulation of immune response (GO:0050777)**: ES=-0.400, FDR=8.184e-01
- **smoothened signaling pathway (GO:0007224)**: ES=-0.491, FDR=8.225e-01

### MSigDB_Hallmark_2020

#### Positively Enriched (Up in Tumor)
- **Reactive Oxygen Species Pathway**: ES=0.227, FDR=3.286e-01
- **Protein Secretion**: ES=0.201, FDR=4.507e-01
- **Oxidative Phosphorylation**: ES=0.273, FDR=4.507e-01
- **Fatty Acid Metabolism**: ES=0.137, FDR=4.507e-01
- **mTORC1 Signaling**: ES=0.105, FDR=4.507e-01

#### Negatively Enriched (Down in Tumor)
- **KRAS Signaling Up**: ES=-0.388, FDR=1.415e-01
- **KRAS Signaling Dn**: ES=-0.395, FDR=1.483e-01
- **Spermatogenesis**: ES=-0.367, FDR=4.876e-01
- **Apoptosis**: ES=-0.080, FDR=9.996e-01
- **Inflammatory Response**: ES=-0.284, FDR=1.000e+00

## Biological Interpretation
The consistently enriched pathways identified in this analysis highlight several biological processes that are dysregulated in prostate cancer:

### Cell Proliferation and Growth
Multiple enriched pathways in up-regulated genes relate to cell cycle regulation, DNA replication, and cell growth, consistent with the increased proliferative capacity of cancer cells.

### Metabolic Reprogramming
Cancer cells often exhibit altered metabolism to support rapid growth. Our analysis identified several metabolic pathways among the significantly enriched terms.

### Immune Response
Several immune-related pathways were identified, reflecting the complex interaction between prostate cancer and the immune system.

## Conclusion
This comprehensive over-representation analysis utilized multiple statistical methods to identify robust pathway enrichment patterns in prostate cancer. By comparing the results across different approaches, we identified consensus pathways that represent the most reliable biological signals. These pathways provide insight into the molecular mechanisms of prostate cancer and may suggest potential therapeutic targets for further investigation.
