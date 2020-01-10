import networkx as nx
import numpy as np
import pandas as pd
import os, pickle, random
import matplotlib.pyplot as plt
from scipy.stats import norm, poisson

VERBOSE = True
try:
    from progress import progress_bar
except:
    VERBOSE = False

GAUSSIAN_NORM = 1e4

class GeneGraph():
    def __init__(self, X:np.array, genes:np.array):
        self.X = np.log(X+1)   # log(Y+1) = k * log(X+1) + \epsilon
        self.weight_mean = np.corrcoef(self.X, rowvar=False) # pearson score

        self.genes = genes
        self.ngene = len(genes)
        self.gene2idx = dict(zip(genes, list(range(self.ngene))))

        # variance of distribution by jackknife
        tmp = np.zeros(shape=[self.ngene, self.ngene, self.ngene])
        for t in range(self.ngene):
            idx = np.array([x for x in range(self.ngene) if x!=t])
            tmp[t, :, :] = np.corrcoef(self.X[idx, :], rowvar=False)
        tmp = tmp - self.weight_mean
        tmp = tmp ** 2
        self.weight_std = np.sum(tmp, axis=0)/(self.ngene)
        self.weight_std = np.maximum(np.absolute(self.weight_std), 0.01)
        # print(self.weight_std)

        self.g = nx.Graph() # network skeleton
        self.g.add_nodes_from(self.genes)

        for idx1, n1 in enumerate(self.genes):
            for idx2, n2 in enumerate(self.genes):
                if idx1 != idx2:
                    if self.weight_std[idx1, idx2] <= 0.01:
                        self.g.add_edge(n1, n2, weight=self.weight_mean[idx1, idx2])

        self.pos = None
        self.approx_cov = self.approx_cov_from_skeleton()

        pass



    def approx_cov_from_skeleton(self):
        # the approximate likelihood matrix generated from network skeleton
        m = np.zeros(shape=[self.ngene, self.ngene])
        for idx1, n1 in enumerate(self.genes):
            for idx2, n2 in enumerate(self.genes):

                # pairwise likelihood   prob = 1 - (1-p1)(1-p2)(1-p3)...
                result = 1
                for a in self.g.neighbors(n1):
                    if a == n2:
                        result *= (1 - self.g[n1][a]['weight'])
                    for b in self.g.neighbors(a):
                        if b not in [n1, a]:
                            if b == n2:
                                result *= (1 - self.g[n1][a]['weight'] * self.g[a][b]['weight'])
                            for c in self.g.neighbors(b):
                                if c not in [n1, a, b]:
                                    if c == n2:
                                        result *= (1 - self.g[n1][a]['weight'] * self.g[a][b]['weight'] * self.g[b][c]['weight'])

                m[idx1][idx2] = 1 - result
        self.approx_cov = m
        return m

    def likelihood(self):
        result = 1.0
        print("Likelihood")
        for idx1, n1 in enumerate(self.genes):
            for idx2, n2 in enumerate(self.genes):
                if idx1 != idx2:
                    prob = norm.pdf(self.approx_cov[idx1, idx2], loc=self.weight_mean[idx1, idx2], scale=self.weight_std[idx1, idx2])
                    result *= max(prob, 1e-4)*1e2
                    # print(result)
        print(result)
        score = result * poisson.pmf(nx.number_of_edges(self.g), POISSON_MU)
        return score

    def position(self, pos=None):
        # assign positions before drawing network
        if pos == None:
            self.pos = nx.spectral_layout(self.g)
        else:
            self.pos = pos

    def print(self, filename="figure/net/tmp.pdf"):

        edges = list(self.g.edges)
        for n1, n2 in edges:
            if self.g[n1][n2]["weight"]<0.3:
                self.g.remove_edge(n1, n2)

        weights = [self.g[u][v]['weight'] for u, v in self.g.edges()]
        plt.figure(figsize=(5, 5))
        nx.draw(self.g, self.pos, with_labels=True, node_size=90, width=weights, node_color="grey")
        plt.savefig(filename)
        plt.show()


N = 1
DELETE_PROB = 0.3
ADD_PROB = 0.3
POISSON_MU = 190


def mcmc(X:np.array, genes:np.array):

    iter_g = GeneGraph(X=X, genes=genes)

    accept = []

    for t in range(N):

        if VERBOSE:
            progress_bar(current=t, total=N,
                         msg='Avg. Accept %d %%' % (np.sum([x[2] for x in accept[-100:]])), epoch=0)
        else:
            if (t+1) % (N/100) ==0:
                print(t+1)

        # Action 1: Delete an edge
        # Action 2: Add an edge
        # Action 2: Change a weight

        rand = random.random()
        accept_rand = random.random()

        if rand <= DELETE_PROB:
            all_edges = list(iter_g.g.edges)
            if len(all_edges) > 0:
                n1, n2 = random.choice(all_edges)
                old_weight = iter_g.g.edges[n1, n2]["weight"]
                mean_weight = iter_g.weight_mean[iter_g.gene2idx[n1], iter_g.gene2idx[n2]]
                std_weight = iter_g.weight_std[iter_g.gene2idx[n1], iter_g.gene2idx[n2]]
                accept_factor = norm.pdf(0, loc=mean_weight, scale=std_weight) / \
                                norm.pdf(old_weight, loc=mean_weight, scale=std_weight)

                # avoid 0
                if poisson.pmf(nx.number_of_edges(iter_g.g), POISSON_MU) < 1e-6:
                    if nx.number_of_edges(iter_g.g) < POISSON_MU:
                        accept_factor *= 0.1
                    else:
                        accept_factor *= 10
                else:
                    accept_factor *= poisson.pmf(nx.number_of_edges(iter_g.g) - 1, POISSON_MU) / \
                                 poisson.pmf(nx.number_of_edges(iter_g.g), POISSON_MU)

                if accept_rand <= accept_factor:
                    iter_g.g.remove_edge(n1, n2)
                    # print(t, "accept delete", accept_factor)
                    accept.append([t, "delete", 1])
                else:
                    accept.append([t, "delete", 0])


        elif rand <= DELETE_PROB + ADD_PROB:
            all_edges = list(nx.non_edges(iter_g.g))
            if len(all_edges) > 0:
                n1, n2 = random.choice(all_edges)

                mean_weight = iter_g.weight_mean[iter_g.gene2idx[n1], iter_g.gene2idx[n2]]
                std_weight = iter_g.weight_std[iter_g.gene2idx[n1], iter_g.gene2idx[n2]]
                new_weight = norm.rvs(loc=mean_weight, scale=1)

                accept_factor = norm.pdf(new_weight, loc=mean_weight, scale=std_weight) / \
                                norm.pdf(0, loc=mean_weight, scale=std_weight)

                # avoid 0
                if poisson.pmf(nx.number_of_edges(iter_g.g), POISSON_MU) < 1e-6:
                    if nx.number_of_edges(iter_g.g) < POISSON_MU:
                        accept_factor *= 10
                    else:
                        accept_factor *= 0.1
                else:
                    accept_factor *= poisson.pmf(nx.number_of_edges(iter_g.g) + 1, POISSON_MU) / \
                                 poisson.pmf(nx.number_of_edges(iter_g.g), POISSON_MU)

                if accept_rand <= accept_factor:
                    iter_g.g.add_edge(n1, n2, weight=new_weight)
                    # print(t, "accept add", accept_factor)
                    accept.append([t, "add", 1])
                else:
                    accept.append([t, "add", 0])

        else:
            all_edges = list(iter_g.g.edges)
            if len(all_edges) > 0:
                n1, n2 = random.choice(all_edges)

                old_weight = iter_g.g.edges[n1, n2]["weight"]
                new_weight = norm.rvs(loc=old_weight, scale=1)

                mean_weight = iter_g.weight_mean[iter_g.gene2idx[n1], iter_g.gene2idx[n2]]
                std_weight = iter_g.weight_std[iter_g.gene2idx[n1], iter_g.gene2idx[n2]]
                accept_factor = norm.pdf(new_weight, loc=mean_weight, scale=std_weight)/\
                                norm.pdf(old_weight, loc=mean_weight, scale=std_weight)

                if accept_rand <= accept_factor:
                    iter_g.g.add_edge(n1, n2, weight=new_weight)
                    # print(t, "accept change", accept_factor)
                    accept.append([t, "change", 1])
                else:
                    accept.append([t, "change", 0])

    print(iter_g.likelihood())

    return iter_g, accept


















if __name__ == "__main__":
    from scvi.dataset import LoomDataset, CsvDataset, Dataset10X, DownloadableAnnDataset
    from scvi.dataset import BrainLargeDataset, CortexDataset, PbmcDataset, RetinaDataset, HematoDataset, CbmcDataset, \
        BrainSmallDataset, SmfishDataset

    save_path = "data/"
    smfish_dataset = SmfishDataset(save_path=save_path)

    g, acc = mcmc(X=smfish_dataset.X, genes=smfish_dataset.gene_names)

    g.print()

    with open("result_{}_{}.pickle".format("osmFISH", str(POISSON_MU)), "wb") as f:
        pickle.dump(obj=[g, acc], file=f)