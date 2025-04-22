Project Overview

This project performs a regression-based analysis to identify transcription factor (TF) and target gene (TG) relationships from gene expression data. The goal is to:

Identify statistically significant relationships between TFs and their potential targets using linear regression.

Compare these predicted relationships with known TF-TG interactions from databases.

Visualize the intersecting known relationships to highlight meaningful biological interactions.

The analysis pipeline involves preprocessing the gene expression data through log2 transformation and Z-score normalization, performing regression for TF-target inference, applying multiple testing correction, and visualizing statistically significant relationships.

Input Files

1. data.csv
Description: Gene expression matrix.
Format: Rows are genes (with gene names as row index), columns are samples.

2. all_tfs.csv
Description: List of transcription factors to test in the regression analysis.
Format: A single-column CSV with header TF

3. known_tf_targets.csv
Description: Known TF-target relationships from curated databases.
Format: Two-column CSV with headers TF and Target.

What the Code Does:

Load and Preprocess Data:
Reads the expression data, TF list, and known relationships.
Applies log2 transformation and Z-score normalization to expression values.

Regression Analysis:
For each TF in the list, performs linear regression against all other genes.
Stores coefficient, p-value, and R-squared values for each TF-target pair.
Applies Benjamini-Hochberg correction for multiple testing.

Filter and Compare with Known Data:
Extracts only those significant TF-target pairs that also exist in the known TF-target database.

Visualization:
Creates a barplot for the top positive and negative known TF-target interactions.
Displays and plots top target genes (positive and negative) for each TF.

Output:
Saves visualizations as PNG files:
known_tf_target_relationships.png
top_pos_neg_targets_per_tf.png
Prints summary of top TF-target relationships to console.

How to Run:
python tf_tg_analysis.py
Ensure all input files are in the correct path or update file paths in the script accordingly.

Requirements:
Python 3.x
pandas
numpy
matplotlib
seaborn
statsmodels
