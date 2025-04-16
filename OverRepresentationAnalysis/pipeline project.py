import os
import warnings
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import gseapy as gp
from scipy.stats import hypergeom, fisher_exact, chi2_contingency, ttest_ind, binomtest
from statsmodels.stats.multitest import multipletests
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score

warnings.filterwarnings('ignore')

# FILE PATHS
DATA_PATH = r"C:\Users\anish\OneDrive\8th SEMESTER\GDA\prostate_exp_cancer_preprocess.csv"
GMT_PATH  = r"C:\Users\anish\OneDrive\8th SEMESTER\GDA\c2.cp.kegg.v7.5.1.symbols.gmt"

##############################################
# Step 1: Data Loading and Preprocessing
##############################################
def load_data(data_path):
    """
    Load dataset and set gene IDs as the index.
    """
    data = pd.read_csv(data_path)
    data = data.set_index(data.columns[0])
    print(f"Loaded data with {data.shape[0]} genes and {data.shape[1]} samples.")
    return data

def filter_low_expression_genes(expr, min_expression=1, min_variance_percentile=10):
    """
    Filter out genes with low mean expression and low variance.
    """
    gene_mean = expr.mean(axis=1)
    gene_var = expr.var(axis=1)
    high_exp_genes = gene_mean[gene_mean > min_expression].index
    var_threshold = np.percentile(gene_var, min_variance_percentile)
    high_var_genes = gene_var[gene_var > var_threshold].index
    filtered_genes = list(set(high_exp_genes) & set(high_var_genes))
    print(f"Filtered genes: from {expr.shape[0]} to {len(filtered_genes)} genes.")
    return expr.loc[filtered_genes]

def identify_highly_expressed_genes(expr, percentile_threshold=90):
    """
    Identify highly expressed genes (top percentile by mean expression).
    """
    gene_means = expr.mean(axis=1)
    threshold = np.percentile(gene_means, percentile_threshold)
    high_expr_genes = gene_means[gene_means >= threshold].index.tolist()
    print(f"Identified {len(high_expr_genes)} highly expressed genes (top {100 - percentile_threshold}%).")
    return high_expr_genes

##############################################
# Additional: Unsupervised Clustering & Visualization
##############################################
def unsupervised_analysis(expr):
    """
    Transpose the expression matrix (samples as rows), standardize,
    run PCA (2 components), and cluster with KMeans.
    Save a PCA scatter plot labeled by clusters.
    """
    expr_T = expr.T
    scaler = StandardScaler()
    scaled_data = scaler.fit_transform(expr_T)
    pca = PCA(n_components=2)
    pca_result = pca.fit_transform(scaled_data)
    pca_df = pd.DataFrame(pca_result, index=expr_T.index, columns=["PC1", "PC2"])
    
    best_score = -1
    best_n = 2
    for n in range(2, 10):
        kmeans = KMeans(n_clusters=n, random_state=42)
        labels = kmeans.fit_predict(pca_df)
        score = silhouette_score(pca_df, labels)
        if score > best_score:
            best_score = score
            best_n = n
    print(f"Optimal clusters: {best_n} (silhouette score: {best_score:.3f})")
    
    kmeans = KMeans(n_clusters=best_n, random_state=42)
    cluster_labels = kmeans.fit_predict(pca_df)
    pca_df['Cluster'] = cluster_labels
    
    plt.figure(figsize=(10,8))
    sns.scatterplot(data=pca_df, x="PC1", y="PC2", hue="Cluster", palette="viridis", s=100)
    plt.title("PCA of Samples with KMeans Clustering")
    plt.xlabel("PC1")
    plt.ylabel("PC2")
    plt.tight_layout()
    plt.savefig("cluster_pca.png", dpi=300, bbox_inches="tight")
    plt.show()
    
    return pca_df, cluster_labels

##############################################
# Extra Visualizations
##############################################
def plot_sample_correlation(expr):
    """
    Plot a heatmap of the sample correlation matrix.
    """
    corr_mat = expr.corr()
    plt.figure(figsize=(10,8))
    sns.heatmap(corr_mat, cmap='vlag', linewidths=0.5)
    plt.title("Sample Correlation Heatmap")
    plt.tight_layout()
    plt.savefig("correlation_heatmap.png", dpi=300, bbox_inches="tight")
    plt.show()

def differential_expression_analysis(expr):
    """
    A simple differential expression analysis by splitting samples into two groups.
    (This is a demo approach if no metadata is provided.)
    Returns a DataFrame with log2 fold change and p-values.
    """
    samples = expr.columns
    mid = len(samples) // 2
    group1 = samples[:mid]
    group2 = samples[mid:]
    results = []
    for gene in expr.index:
        vals1 = expr.loc[gene, group1]
        vals2 = expr.loc[gene, group2]
        if vals1.mean() <= 0 or vals2.mean() <= 0:
            continue
        log2fc = np.log2(vals1.mean() / vals2.mean())
        t_stat, p_val = ttest_ind(vals1, vals2, equal_var=False)
        results.append({'gene': gene, 'log2fc': log2fc, 'p_value': p_val})
    de_df = pd.DataFrame(results)
    de_df['adjusted_p_value'] = multipletests(de_df['p_value'], method='fdr_bh')[1]
    return de_df

def plot_volcano(de_df, fc_col='log2fc', padj_col='adjusted_p_value'):
    """
    Plot a volcano plot from differential expression results.
    """
    df = de_df.copy()
    df['-log10(padj)'] = -np.log10(df[padj_col] + 1e-10)
    plt.figure(figsize=(10,8))
    sns.scatterplot(data=df, x=fc_col, y='-log10(padj)', alpha=0.7, edgecolor=None)
    plt.axhline(-np.log10(0.05), color='red', linestyle='--', label='adj p = 0.05')
    plt.xlabel("Log2 Fold Change")
    plt.ylabel("-Log10 Adjusted p-value")
    plt.title("Volcano Plot of Differential Expression")
    plt.legend()
    plt.tight_layout()
    plt.savefig("volcano_plot.png", dpi=300, bbox_inches="tight")
    plt.show()

##############################################
# Step 2: Load GMT File with Gene Sets
##############################################
def load_gmt(gmt_path):
    """
    Load gene sets from the GMT file.
    """
    pathways = {}
    with open(gmt_path, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            term = parts[0]
            desc = parts[1]
            genes = set(parts[2:])
            pathways[term] = {'description': desc, 'genes': genes}
    print(f"Loaded {len(pathways)} pathways from GMT file.")
    return pathways

##############################################
# Step 3: Define ORA Test Functions
##############################################
def hypergeometric_test(gene_list, pathway_genes, background):
    gene_set = set(gene_list)
    pathway_set = set(pathway_genes)
    background_set = set(background)
    M = len(background_set)
    n = len(pathway_set & background_set)
    N = len(gene_set & background_set)
    k = len(gene_set & pathway_set)
    if k == 0 or n == 0 or N == 0:
        return 1.0, 0, []
    p_val = hypergeom.sf(k - 1, M, n, N)
    fold_enrichment = (k/N) / (n/M) if (n > 0 and N > 0) else 0
    return p_val, fold_enrichment, list(gene_set & pathway_set)

def fishers_exact_test(gene_list, pathway_genes, background):
    gene_set = set(gene_list)
    pathway_set = set(pathway_genes)
    background_set = set(background)
    a = len(gene_set & pathway_set)
    b = len(pathway_set & background_set) - a
    c = len(gene_set) - a
    d = len(background_set) - (a + b + c)
    if a == 0:
        return 1.0, 0, []
    table = [[a, b], [c, d]]
    _, p_val = fisher_exact(table, alternative='greater')
    fold_enrichment = (a / len(gene_set)) / ((a+b) / len(background_set)) if len(gene_set)>0 and (a+b)>0 else 0
    return p_val, fold_enrichment, list(gene_set & pathway_set)

def binomial_test_method(gene_list, pathway_genes, background):
    gene_set = set(gene_list)
    pathway_set = set(pathway_genes)
    background_set = set(background)
    overlap = len(gene_set & pathway_set)
    if overlap == 0:
        return 1.0, 0, []
    p_success = len(pathway_set & background_set) / len(background_set)
    result = binomtest(overlap, n=len(gene_set), p=p_success, alternative='greater')
    p_val = result.pvalue
    fold_enrichment = (overlap / len(gene_set)) / p_success if p_success > 0 else 0
    return p_val, fold_enrichment, list(gene_set & pathway_set)

def chi_squared_test_method(gene_list, pathway_genes, background):
    gene_set = set(gene_list)
    pathway_set = set(pathway_genes)
    background_set = set(background)
    a = len(gene_set & pathway_set)
    b = len(pathway_set & background_set) - a
    c = len(gene_set) - a
    d = len(background_set) - (a + b + c)
    if a == 0:
        return 1.0, 0, []
    table = [[a, b], [c, d]]
    _, p_val, _, _ = chi2_contingency(table)
    fold_enrichment = (a / len(gene_set)) / ((a+b) / len(background_set)) if len(gene_set)>0 and (a+b)>0 else 0
    return p_val, fold_enrichment, list(gene_set & pathway_set)

##############################################
# Step 4: Run All ORA Methods on Gene List
##############################################
def run_all_ora_methods(gene_list, gmt_path, background):
    pathways = load_gmt(gmt_path)
    hyper_results, fisher_results = [], []
    binomial_results, chi2_results = [], []
    
    for term, info in pathways.items():
        genes_in_pathway = info['genes']
        p_hyper, fe_hyper, overlap_hyper = hypergeometric_test(gene_list, genes_in_pathway, background)
        hyper_results.append({
            'term_id': term,
            'term_name': info['description'],
            'p_value': p_hyper,
            'fold_enrichment': fe_hyper,
            'overlap_count': len(overlap_hyper),
            'overlap_genes': ','.join(overlap_hyper)
        })
        p_fisher, fe_fisher, overlap_fisher = fishers_exact_test(gene_list, genes_in_pathway, background)
        fisher_results.append({
            'term_id': term,
            'term_name': info['description'],
            'p_value': p_fisher,
            'fold_enrichment': fe_fisher,
            'overlap_count': len(overlap_fisher),
            'overlap_genes': ','.join(overlap_fisher)
        })
        p_binom, fe_binom, overlap_binom = binomial_test_method(gene_list, genes_in_pathway, background)
        binomial_results.append({
            'term_id': term,
            'term_name': info['description'],
            'p_value': p_binom,
            'fold_enrichment': fe_binom,
            'overlap_count': len(overlap_binom),
            'overlap_genes': ','.join(overlap_binom)
        })
        p_chi2, fe_chi2, overlap_chi2 = chi_squared_test_method(gene_list, genes_in_pathway, background)
        chi2_results.append({
            'term_id': term,
            'term_name': info['description'],
            'p_value': p_chi2,
            'fold_enrichment': fe_chi2,
            'overlap_count': len(overlap_chi2),
            'overlap_genes': ','.join(overlap_chi2)
        })
    
    hyper_df = pd.DataFrame(hyper_results)
    fisher_df = pd.DataFrame(fisher_results)
    binomial_df = pd.DataFrame(binomial_results)
    chi2_df = pd.DataFrame(chi2_results)
    
    for df in [hyper_df, fisher_df, binomial_df, chi2_df]:
        df['adjusted_p_value'] = multipletests(df['p_value'], method='fdr_bh')[1]
        df.sort_values('adjusted_p_value', inplace=True)
    
    print("ORA methods completed.")
    return {
        'Hypergeometric': hyper_df,
        'Fisher': fisher_df,
        'Binomial': binomial_df,
        'ChiSquared': chi2_df
    }

##############################################
# Step 5: Generate Summary Report
##############################################
def generate_summary_report(ora_results, report_file="ora_summary_report.md"):
    report_lines = []
    report_lines.append("# Summary Report: Over-Representation Analysis\n")
    report_lines.append("## Overview\n")
    report_lines.append("This report summarizes ORA results using four tests:\n- Hypergeometric\n- Fisher's Exact\n- Binomial\n- Chi-Squared\n")
    
    for method, df in ora_results.items():
        enriched = df[df['adjusted_p_value'] < 0.05]
        report_lines.append(f"### {method} Test")
        report_lines.append(f"- Total pathways tested: {df.shape[0]}")
        report_lines.append(f"- Enriched pathways (adj p < 0.05): {enriched.shape[0]}\n")
        report_lines.append("#### Top 5 Enriched Pathways:")
        top5 = enriched.head(5)
        for idx, row in top5.iterrows():
            report_lines.append(f"- **{row['term_id']}** ({row['term_name']}): adj p = {row['adjusted_p_value']:.2e}, fold enrichment = {row['fold_enrichment']:.2f}, overlap = {row['overlap_count']}")
        report_lines.append("\n")
        
    with open(report_file, 'w') as f:
        f.write("\n".join(report_lines))
    print(f"Summary report generated: {report_file}")

##############################################
# Main Execution: Run Pipeline and Visualize Everything
##############################################
if __name__ == "__main__":
    # Step 1: Load data and filter genes
    data = load_data(DATA_PATH)
    filtered_data = filter_low_expression_genes(data, min_expression=1, min_variance_percentile=10)
    
    # Unsupervised clustering and PCA visualization
    pca_df, cluster_labels = unsupervised_analysis(filtered_data)
    
    # Plot sample correlation heatmap
    plot_sample_correlation(filtered_data)
    
    # Optional: Differential Expression Analysis to create a Volcano Plot
    de_results = differential_expression_analysis(filtered_data)
    if not de_results.empty:
        plot_volcano(de_results)
    
    # Step 2: Identify highly expressed genes (to serve as input gene list for ORA)
    gene_list = identify_highly_expressed_genes(filtered_data, percentile_threshold=90)
    background = filtered_data.index.tolist()
    
    # Step 3: Run ORA methods on the selected gene list
    ora_results = run_all_ora_methods(gene_list, GMT_PATH, background)
    for method, df in ora_results.items():
        filename = f"{method}_ora_results.csv"
        df.to_csv(filename, index=False)
        print(f"Saved {method} results to {filename}")
    
    # Step 4: Generate summary report
    generate_summary_report(ora_results, report_file="ora_summary_report.md")
    
    print("All analyses and visualizations completed successfully!")
