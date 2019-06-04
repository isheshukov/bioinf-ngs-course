import fire
import networkx as nx
from Bio import SeqIO
from itertools import product
from more_itertools import windowed
import dataclasses
import os 
from tqdm import tqdm
from collections import Counter
from numpy import mean
from copy import deepcopy

def graph2fasta(G, name):
    with open(name, "w") as file:
        for e in G.edges(data=True):
            file.write(f"> {e[0]} -> {e[1]}, length = {e[2]['length']}, mean_cov = {e[2]['mean_cov']}\n")
            file.write(f"{e[2]['seq']}\n")

def seq2kmers(seq, k):
    return map(lambda x: ''.join(x) if None not in x else "", windowed(seq, k))

def fasta2graph(fasta, type, k):
    coverage = Counter() 
    G = nx.DiGraph()

    seqs = []
    with open(fasta, "r") as handle:
        for record in tqdm(SeqIO.parse(handle, type), desc="Reading sequences"):
            seqs.append(record.seq)
            seqs.append(record.reverse_complement().seq)

    # Adding kmers as nodes
    for seq in tqdm(seqs, desc="Creating nodes"):
        for kmer in seq2kmers(seq, k):
            if kmer != "":
              G.add_node(kmer)
              coverage[kmer] += 1
    
    # Adding k+1-mers as edges
    for seq in tqdm(seqs, desc="Creating edges"):
        for kmer in seq2kmers(seq, k+1):
            if kmer != "":
                l = kmer[:-1]
                r = kmer[1:]
                if (l in G) and (r in G):
                    G.add_edge(l, r, seq=kmer, length=len(kmer))

    #for l, r in tqdm(product(G, repeat=2), desc="Creating edges", total=len(G)**2):
    #    if l[1:] == r[:-1]:
    #        G.add_edge(l, r, seq=l+r[-1], length=len(l+r[-1]))


    return G, coverage

def compress_graph(G, k):
    nodes_to_remove = []

    for node in tqdm(G, desc="Compressing"):
        if G.in_degree(node) == G.out_degree(node) == 1:
            in_node, _, in_data = list(G.in_edges(node, data=True))[0]
            _, out_node, out_data = list(G.out_edges(node, data=True))[0]

            new_seq = in_data['seq'] + out_data['seq'][k:]
            G.add_edge(in_node, out_node, seq=new_seq, length=len(new_seq))

            G.remove_edge(in_node, node)
            G.remove_edge(node, out_node)
            nodes_to_remove.append(node)


    G.remove_nodes_from(nodes_to_remove)

def compute_mean_coverage(G, coverage, k):
    def seq_to_mean(seq, k):
        m = c = 0
        for i in seq2kmers(seq, k):
            m += coverage[i]
            c += 1

        return m/c

    labels_dict = {}
    mean_cov = {}
    #length = {}
    for e in tqdm(G.edges(data=True), desc="Computing mean coverage"):
        mean_cov[(e[0],e[1])] = seq_to_mean(e[2]['seq'], k)

        labels_dict[(e[0],e[1])] = f"mean_cov: {float(seq_to_mean(e[2]['seq'], k)):0.3}\n length: {e[2]['length']}"
    
    nx.set_edge_attributes(G, labels_dict, 'label')
    nx.set_edge_attributes(G, mean_cov, 'mean_cov')
    #nx.set_edge_attributes(G, length, 'length')

def remove_bad_edges(G, coverage, k):
    len_threshold = 2*k
    cov_threshold = 33

    edges_to_remove = []
    for e in tqdm(G.edges(data=True), desc="Finding bad edges"):
        if e[2]['mean_cov'] < cov_threshold or e[2]['length'] < len_threshold:
            edges_to_remove.append(e) 

    for i in edges_to_remove:
            G.remove_edge(i[0], i[1])

    nodes_to_remove = []
    for node in tqdm(G, desc="Removing orphaned nodes"):
        if G.in_degree(node) == G.out_degree(node) == 0:
            nodes_to_remove.append(node)

    G.remove_nodes_from(nodes_to_remove)

def remove_dead_ends(G, coverage, k):
    len_threshold = 2*k 
    cov_threshold = 33 
    edges_to_remove = []
    for node in tqdm(G, desc="Compressing"):
        if G.out_degree(node) == 0:
            in_edges = list(G.in_edges(node, data=True))
            if not in_edges:
                break

            in_node, _, in_data = in_edges[0]

            edge_data = G.get_edge_data(in_node, node)
            if edge_data['mean_cov'] < cov_threshold or edge_data['length'] < len_threshold:
                edges_to_remove.append((in_node, node))

    G.remove_edges_from(edges_to_remove)

    nodes_to_remove = []
    for node in tqdm(G, desc="Removing orphaned nodes"):
        if G.in_degree(node) == G.out_degree(node) == 0:
            nodes_to_remove.append(node)

    G.remove_nodes_from(nodes_to_remove)

def main(filename, k, type="fasta"):
    """Runs de Brujin graph building algorithm """

    basename = os.path.splitext(filename)[0]
    extension = os.path.splitext(filename)[1][1:]

    type = type if extension not in ['fasta', 'fastq'] else extension
    dbg, coverage = fasta2graph(filename, type, k)
    compress_graph(dbg, k)

    dbg = nx.relabel.convert_node_labels_to_integers(dbg)
    compute_mean_coverage(dbg, coverage, k)
    dbg_bad = deepcopy(dbg)
    dbg_end = deepcopy(dbg)


    graph2fasta(dbg, basename + ".graph.fasta")
    nx.set_edge_attributes(dbg, '', 'seq')
    nx.nx_agraph.write_dot(dbg, basename + ".dot")

    remove_bad_edges(dbg_bad, coverage, k)
    nx.set_edge_attributes(dbg_bad, '', 'seq')
    graph2fasta(dbg_bad, basename + ".no_bad_edges.graph.fasta")
    nx.nx_agraph.write_dot(dbg_bad, basename + ".no_bad_edges.dot")

    remove_dead_ends(dbg_end, coverage, k)
    nx.set_edge_attributes(dbg_end, '', 'seq')
    graph2fasta(dbg_end, basename + ".no_dead_ends.graph.fasta")
    nx.nx_agraph.write_dot(dbg_end, basename + ".no_dead_ends.dot")

    # We need to strip seq label, otherwise dot will fail


    return

if __name__ == '__main__':
  fire.Fire(main)
