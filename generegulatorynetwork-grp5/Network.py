import pandas as pd
import numpy as np
from sklearn.linear_model import LinearRegression
from sklearn.linear_model import Lasso
from sklearn.preprocessing import StandardScaler
import networkx as nx
import matplotlib.pyplot as plt

class Gene_Network:



    def __init__(self,csv_file,tfs_file):
        self.data = pd.read_csv(csv_file)
        self.df = pd.read_csv(tfs_file)
        self.x = pd.read_csv('try.csv', index_col = 0)
        x = self.x
        scaler = StandardScaler()
        x_scaled = pd.DataFrame(scaler.fit_transform(x.T), index=x.columns, columns=x.index).T
        col_0 = self.data.columns[0]
        self.data.rename(columns={col_0: "Gene"}, inplace=True)
        tfs = self.df['Transcription_factors']
        t = tfs.to_numpy()
        proteins = self.data['Gene']
        p = proteins.to_numpy()
        self.match = []
        for i in range(len(t)):
            tf = t[i]
            for j in range(len(proteins)):
                if tf == proteins[j]:
                    self.match.append(tf)
        self.tf_data = x_scaled.loc[self.match]



        
    def network_one(self,threshold_p=0.7,threshold_n=-0.7):
        x = self.x
        tf_data = self.x.loc[self.match]
        self.x.drop(self.match, axis = 0, inplace = True)
        scaler = StandardScaler()
        tf_data_scaled = pd.DataFrame(scaler.fit_transform(tf_data.T), index=tf_data.columns, columns=tf_data.index).T
        x_scaled = pd.DataFrame(scaler.fit_transform(x.T), index=x.columns, columns=x.index).T
        corr_genes_pos = []
        corr_genes_neg = []
        for i in range(len(self.match)):
            positive_genes = [self.match[i]]
            negative_genes = [self.match[i]]
            print(self.match[i])
            correlation_matrix = x_scaled.T.corrwith(tf_data_scaled.iloc[i], axis=0)
            positive_genes.extend(gene for gene, value in correlation_matrix.items() if value > threshold_p)
            negative_genes.extend(gene for gene, value in correlation_matrix.items() if value < threshold_n)
            corr_genes_pos.append(positive_genes)
            corr_genes_neg.append(negative_genes)
        self.print_clean(corr_genes_pos,corr_genes_neg)



        
    def network_two(self,threshold_p = 0.7, threshold_n= -0.1):
        model = LinearRegression()
        x_t = self.x.T
        tf_data_t = self.tf_data.T
        gene_names = x_t.columns
        tf_names = tf_data_t.columns
        x_train = tf_data_t
        w_of_tfs = []
        i = 0;
        for gene_x in gene_names:
            y_train = x_t[gene_x]
            model.fit(x_train,y_train)
            w = model.coef_
            threshold = 0.5
            w_of_tfs.append(w)
        linear_genes_pos = []
        linear_genes_neg = []
        match = self.match
        for i in range(len(w_of_tfs[0])):
            positive_genes = [match[i]]
            negative_genes = [match[i]]
            for j in range(len(w_of_tfs)):
                if w_of_tfs[j][i]>threshold_p:
                    positive_genes.append(gene_names[j])
                if w_of_tfs[j][i]<threshold_n:
                    negative_genes.append(gene_names[j])
            linear_genes_pos.append(positive_genes)
            linear_genes_neg.append(negative_genes)
        self.print_clean(linear_genes_pos,linear_genes_neg)




    
    def network_three(self,threshold_p=0.7,threshold_n=-0.7):
        lasso = Lasso(alpha=0.1)
        x_t = self.x.T
        tf_data_t = self.tf_data.T
        gene_names = x_t.columns
        tf_names = tf_data_t.columns
        x_train = tf_data_t
        w_of_tfs = []
        i = 0;
        match = self.match
        for gene_x in gene_names:
            y_train = x_t[gene_x]
            lasso.fit(x_train,y_train)
            w = lasso.coef_
            threshold = 0.5
            w_of_tfs.append(w)
        lasso_genes_pos = []
        lasso_genes_neg = []
        for i in range(len(w_of_tfs[0])):
            positive_genes = [match[i]]
            negative_genes = [match[i]]
            for j in range(len(w_of_tfs)):
                if w_of_tfs[j][i]>threshold_p:
                    positive_genes.append(gene_names[j])
                if w_of_tfs[j][i]<threshold_n:
                    negative_genes.append(gene_names[j])
            lasso_genes_pos.append(positive_genes)
            lasso_genes_neg.append(negative_genes)
        self.print_clean(lasso_genes_pos,lasso_genes_neg)
        return lasso_genes_pos, lasso_genes_neg

        
    def print_clean(self,pos,neg):
        print('--------------------------------------------------------------Positive Relation ----------------------------------------------------------------------   ')
        for i in range(len(self.match)):
            print(f'genes positively related with {pos[i][0]} transcription factors are')
            for j in pos[i]:
                if j==pos[i][0]:
                    continue
                print(f'{j}  ',end = '')
            print()
            print()
        print('--------------------------------------------------------------Negative Relation ----------------------------------------------------------------------   ')
        for i in range(len(self.match)):
            print(f'genes negatively related with {neg[i][0]} transcription factors are')
            for j in neg[i]:
                if j== neg[i][0]:
                    continue
                print(f'{j}  ',end = '')
            print()
            print()
        a = pos
        b = neg
        G = nx.DiGraph()

# Add all TFs first
        tfs = set([group[0] for group in a] + [group[0] for group in b])
        for tf in tfs:
            G.add_node(tf, role='TF')

# Add edges from a
        for group in a:
            src = group[0]
            for target in group[1:]:
                G.add_node(target, role='gene')
                G.add_edge(src, target, source='a')

# Add edges from b
        for group in b:
            src = group[0]
            for target in group[1:]:
                G.add_node(target, role='gene')
                if not G.has_edge(src, target):
                    G.add_edge(src, target, source='b')

# Separate TFs and genes for shell layout
        tf_nodes = [n for n, d in G.nodes(data=True) if d['role'] == 'TF']
        gene_nodes = [n for n, d in G.nodes(data=True) if d['role'] == 'gene']
        pos = nx.shell_layout(G, [tf_nodes, gene_nodes])

# Node colors
        node_colors = ['lightblue' if G.nodes[n]['role'] == 'TF' else 'lightcoral' for n in G.nodes()]

# Edge colors
        edge_colors = []
        for u, v in G.edges():
            edge_colors.append('green' if G[u][v].get('source') == 'a' else 'purple')

# Draw clearer plot
        plt.figure(figsize=(14, 10))
        nx.draw_networkx_nodes(G, pos, node_color=node_colors, node_size=1000, edgecolors='black')
        nx.draw_networkx_labels(G, pos, font_size=9)
        nx.draw_networkx_edges(G, pos, edge_color=edge_colors, arrows=True, arrowstyle='->', width=2)

# Optional: Legend
        import matplotlib.patches as mpatches
        green_patch = mpatches.Patch(color='green', label='Positive regulation')
        purple_patch = mpatches.Patch(color='purple', label='Negative regulation')
        blue_patch = mpatches.Patch(color='lightblue', label='Transcription Factor')
        pink_patch = mpatches.Patch(color='lightcoral', label='Target Gene')
        plt.legend(handles=[green_patch, purple_patch, blue_patch, pink_patch], loc='upper left')

        plt.title("Clean Gene Regulatory Network", fontsize=14)
        plt.axis('off')
        plt.tight_layout()
        plt.show()





