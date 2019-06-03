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

def graph2fasta(G, name):
    with open(name, "w") as file:
        for e in G.edges(data=True):
            file.write(f"> {e[0]} -> {e[1]}\n")
            file.write(f"{e[2]['seq']}\n")

def seq2kmers(seq, k):
    return map(lambda x: ''.join(x) if None not in x else "", windowed(seq, k))

def fasta2graph(fasta, type, k):
    coverage = Counter() 
    G = nx.MultiDiGraph()

    seqs = []
    with open(fasta, "r") as handle:
        for record in tqdm(SeqIO.parse(handle, type), desc="Reading sequences"):
            seqs.append(record.seq)
            seqs.append(record.reverse_complement().seq)

    for seq in tqdm(seqs, desc="Creating nodes"):
        for kmer in seq2kmers(seq, k):
            if kmer != "":
              G.add_node(kmer)
              coverage[kmer] += 1

    for l, r in tqdm(product(G, repeat=2), desc="Creating edges", total=len(G)**2):
        if l[1:] == r[:-1]:
            G.add_edge(l, r, seq=l+r[-1], length=len(l+r[-1]))


    return G, coverage

def compress_graph(G, k):
    nodes_to_remove = []

    for node in tqdm(G, desc="Compressing"):
        if G.in_degree(node) == G.out_degree(node) == 1:
            in_node, _, in_data = list(G.in_edges(node, data=True))[0]
            _, out_node, out_data = list(G.out_edges(node, data=True))[0]

            G.remove_edge(in_node, node)
            G.remove_edge(node, out_node)
            nodes_to_remove.append(node)

            new_seq = in_data['seq'] + out_data['seq']
            print(len(new_seq))
            print('-----')
            G.add_edge(in_node, out_node, seq=new_seq, length=len(new_seq))

    G.remove_nodes_from(nodes_to_remove)

def compute_mean_coverage(G, coverage, k):
    def seq_to_mean(seq, k):
        m = 0
        c = 0
        for i in seq2kmers(seq, k):
            m += coverage[i]
            c += 1

        return m/c

    labels_dict = {}
    mean_cov = {}
    #length = {}

    for e in tqdm(G.edges(data=True, keys=True), desc="Computing mean coverage"):
        mean_cov[(e[0],e[1],e[2])] = seq_to_mean(e[3]['seq'], k)
        #length[(e[0],e[1],e[2])] = len(e[3]['seq'])

        labels_dict[(e[0],e[1],e[2])] = f"mean_cov: {seq_to_mean(e[3]['seq'], k):0.3}\n length: {e[3]['length']}"
    
    nx.set_edge_attributes(G, labels_dict, 'label')
    nx.set_edge_attributes(G, mean_cov, 'mean_cov')
    #nx.set_edge_attributes(G, length, 'length')

def remove_bad_edges(G, coverage, k):
    remove_threshold = 2*k

    edges_to_remove = []
    for e in tqdm(G.edges(data=True, keys=True), desc="Finding bad edges"):
        #if e[3]['mean_cov'] < 4 or e[3]['length'] < remove_threshold:
        if e[3]['length'] < remove_threshold:
            edges_to_remove.append(e) 

    G.remove_edges_from(edges_to_remove)

    nodes_to_remove = []
    for node in tqdm(G, desc="Removing orphaned nodes"):
        if G.in_degree(node) == G.out_degree(node) == 0:
            nodes_to_remove.append(node)

    G.remove_nodes_from(nodes_to_remove)

def remove_dead_ends(G, coverage, k):
    threshold = 2*k
    edges_to_remove = []
    for node in tqdm(G, desc="Compressing"):
        if G.out_degree(node) == 0:
            in_node, _, in_data = list(G.in_edges(node, data=True))[0]

            edge_data = G.get_edge_data(in_node, node)
            #print(edge_data)
            if 'length' in edge_data:
                if (edge_data['length'] < threshold):
                    edges_to_remove.append((in_node, node))

    G.remove_edges_from(edges_to_remove)

def main(filename, k, type="fasta"):
    """Runs de Brujin graph building algorithm """

    basename = os.path.splitext(filename)[0]
    extension = os.path.splitext(filename)[1][1:]

    type = type if extension not in ['fasta', 'fastq'] else extension
    dbg, coverage = fasta2graph(filename, type, k)
    compress_graph(dbg, k)

    compute_mean_coverage(dbg, coverage, k)
    remove_bad_edges(dbg, coverage, k)
    #remove_dead_ends(dbg, coverage, k)

    dbg = nx.relabel.convert_node_labels_to_integers(dbg)
    #dbg = compress_graph(dbg)

    graph2fasta(dbg, basename + ".graph.fasta")

    # We need to strip seq label, otherwise dot will fail

    nx.set_edge_attributes(dbg, [''], 'seq')
    nx.nx_agraph.write_dot(dbg, basename + ".dot")

    return

if __name__ == '__main__':
  fire.Fire(main)