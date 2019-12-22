import networkx as nx
import numpy as np
import pandas as pd
import os, pickle
import matplotlib.pyplot as plt



class GeneGraph():
    def __init__(self, X:np.array, genes:np.array, cells):
        self.X = X
        self.cov = np.cov(X.T)
        self.genes = genes
        self.cells = cells
        self.g = nx.Graph()
        self.g.add_nodes_from(self.genes)
        self.ngene = len(genes)
        self.pos = None
        pass

    def pair_likelihood(self, n1, n2, G=None):
        if G == None:
            G = self.g
        # prob = 1 - (1-p1)(1-p2)(1-p3)...
        result = 1
        for a in G.neighbors(n1):
            if a == n2:
                result *= (1-G[n1][a]['weight'])
            for b in G.neighbors(a):
                if b not in [n1, a]:
                    if b == n2:
                        result *= (1-G[n1][a]['weight']*G[a][b]['weight'])
                    for c in G.neighbors(b):
                        if c not in [n1, a, b]:
                            if c == n2:
                                result *= (1-G[n1][a]['weight']*G[a][b]['weight']*G[b][c]['weight'])
        return 1-result

    def approx_cov(self):
        m = np.zeros(shape=[self.ngene, self.ngene])
        for n1 in range(self.ngene):
            for n2 in range(n1, self.ngene):
                m[n1][n2] = self.pair_likelihood(self.genes[n1], self.genes[n2])
        return m

    def likelihood(self, G):
        
        pass

    def position(self, pos=None):
        if pos == None:
            self.pos = nx.fruchterman_reingold_layout(self.g)
        else:
            self.pos = pos


    def print(self):
        pos = nx.fruchterman_reingold_layout(self.g)

        weights = [self.g[u][v]['weight'] / 100 for u, v in self.g.edges()]
        plt.figure(figsize=(20, 20))
        nx.draw(self.g, pos, with_labels=True, node_size=90, width=weights, node_color="grey")
        plt.savefig("figure/net/tmp.pdf")
        plt.show()

















if __name__ == "__main__":
    from scvi.dataset import LoomDataset, CsvDataset, Dataset10X, DownloadableAnnDataset
    from scvi.dataset import BrainLargeDataset, CortexDataset, PbmcDataset, RetinaDataset, HematoDataset, CbmcDataset, \
        BrainSmallDataset, SmfishDataset

    save_path = "data/"
    smfish_dataset = SmfishDataset(save_path=save_path)
    gene_graph = GeneGraph(X=smfish_dataset.X, genes=smfish_dataset.gene_names, cells=smfish_dataset.labels)
    gene_graph.print()

    print(gene_graph.approx_cov())