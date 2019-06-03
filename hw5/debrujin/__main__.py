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
    return map(lambda x: ''.join(x), windowed(seq, k))

def fasta2graph(fasta, k):
    coverage = Counter() 
    G = nx.MultiDiGraph()

    seqs = []
    with open(fasta, "r") as handle:
        for record in tqdm(SeqIO.parse(handle, "fasta"), desc="Reading sequences"):
            seqs.append(record.seq)
            seqs.append(record.reverse_complement().seq)

    for seq in tqdm(seqs, desc="Creating nodes"):
        for kmer in seq2kmers(seq, k):
            G.add_node(kmer[:-1])
            G.add_node(kmer[1:])
            coverage[kmer] += 1

    for l, r in tqdm(product(G, repeat=2), desc="Creating edges"):
        if l[1:] == r[:-1]:
            G.add_edge(l, r, seq=l+r[-1])


    return G, coverage

def compress_graph(G, k):
    nodes_to_remove = []
    num_of_nodes = len(G)

    while True:
        for node in tqdm(G, desc="Compressing"):
            if G.in_degree(node) == G.out_degree(node) == 1:
                in_node, _, in_data = list(G.in_edges(node, data=True))[0]
                _, out_node, out_data = list(G.out_edges(node, data=True))[0]

                G.remove_edge(in_node, node)
                G.remove_edge(node, out_node)
                nodes_to_remove.append(node)
                G.add_edge(in_node, out_node, seq=in_data['seq']+out_data['seq'])

        G.remove_nodes_from(nodes_to_remove)
        if num_of_nodes == len(G):
            break
        else:
            num_of_nodes = len(G)

def compute_mean_coverage(G, coverage, k):
    def seq_to_mean(seq, k):
        return mean(list(map(lambda x: coverage[''.join(x)], windowed(seq, k))))

    mean_dict = {}

    for e in tqdm(G.edges(data=True, keys=True), desc="Computing mean coverage"):
        mean_dict[(e[0],e[1],e[2])] = {'mean cov': seq_to_mean(e[3]['seq'], k)}
    
    nx.set_edge_attributes(G, mean_dict, 'mean cov')

def main(fasta, k=3):
    """Runs de Brujin graph building algorithm """

    basename = os.path.splitext(fasta)[0]

    dbg, coverage = fasta2graph(fasta, k)
    compress_graph(dbg, k)
    compute_mean_coverage(dbg, coverage, k)

    dbg = nx.relabel.convert_node_labels_to_integers(dbg)
    #dbg = compress_graph(dbg)

    graph2fasta(dbg, basename + ".graph.fasta")

    # We need to strip seq label, otherwise dot will file

    nx.set_edge_attributes(dbg, [''], 'seq')
    nx.nx_agraph.write_dot(dbg, basename + ".dot")

    return

if __name__ == '__main__':
  fire.Fire(main)