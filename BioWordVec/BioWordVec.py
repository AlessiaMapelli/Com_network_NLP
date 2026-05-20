# coding=utf-8
'''
Reference implementation for joint learning word vectors from biomedical text and MeSH graph,
using FastText and modern Node2Vec (PyPI version).
'''

import subprocess
import sys

# Safe imports after ensuring packages are installed
import argparse
import networkx as nx
from node2vec import Node2Vec
from gensim.models import FastText
import gzip
import pickle as pkl


def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument('--input_corpus', default='./data/pubmed_sample',
                        help='Input biomedical corpus (plain text, one sentence per line)')

    parser.add_argument('--input_mesh', default='./data/MeSH_graph.edgelist',
                        help='Input MeSH knowledge graph in edgelist format')

    parser.add_argument('--input_dic', default='./data/MeSH_dic.pkl.gz',
                        help='Pickled MeSH dictionary (gzip)')

    parser.add_argument('--output_model', default='./pubmed_mesh_test',
                        help='Output FastText model path')

    parser.add_argument('--output_bin', default='./pubmed_mesh_test.bin',
                        help='Output word vector binary file')

    parser.add_argument('--dimensions', type=int, default=200,
                        help='Embedding size')

    parser.add_argument('--walk_length', type=int, default=50,
                        help='Length of each random walk')

    parser.add_argument('--num_walks', type=int, default=10,
                        help='Number of walks per node')

    parser.add_argument('--window_size', type=int, default=5,
                        help='Context window size')

    parser.add_argument('--epochs', type=int, default=5,
                        help='Number of training epochs')

    parser.add_argument('--min_count', type=int, default=5,
                        help='Minimum count for words in FastText')

    parser.add_argument('--sg', type=int, default=1,
                        help='Use skip-gram (1) or CBOW (0)')

    parser.add_argument('--workers', type=int, default=4,
                        help='Number of worker threads')

    parser.add_argument('--p', type=float, default=1.0,
                        help='Return hyperparameter for node2vec')

    parser.add_argument('--q', type=float, default=1.0,
                        help='Inout hyperparameter for node2vec')

    return parser.parse_args()


class JointCorpus:
    def __init__(self, mesh_walks, pubmed_path):
        self.mesh_walks = mesh_walks
        self.pubmed_path = pubmed_path

    def __iter__(self):
        for walk in self.mesh_walks:
            yield walk
        with open(self.pubmed_path, 'r') as f:
            for line in f:
                yield line.strip().split()


def main(args):
    # Load MeSH dictionary
    with gzip.open(args.input_dic, 'rb') as f:
        mesh_dict = pkl.load(f)

    # Load MeSH graph
    print("📥 Reading graph...")
    G = nx.read_edgelist(args.input_mesh, nodetype=str, create_using=nx.Graph())

    # Initialize and run node2vec walks
    print("🚶 Generating node2vec walks...")
    node2vec = Node2Vec(
        G,
        dimensions=args.dimensions,
        walk_length=args.walk_length,
        num_walks=args.num_walks,
        workers=args.workers,
        p=args.p,
        q=args.q,
        quiet=True
    )
    walks = node2vec.walks  # List[List[str]]

    # Convert nodes in walks using mesh_dict
    print("🧠 Mapping walks to MeSH terms...")
    mapped_walks = []
    for walk in walks:
        mapped_walks.append([mesh_dict[node] for node in walk if node in mesh_dict])

    # Combine walks and corpus for joint training
    print("⚙️ Training FastText model...")
    model = FastText(
        sentences=JointCorpus(mapped_walks, args.input_corpus),
        vector_size=args.dimensions,
        window=args.window_size,
        min_count=args.min_count,
        workers=args.workers,
        sg=args.sg,
        epochs=args.epochs
    )

    # Save outputs
    model.save(args.output_model)
    print(f"✅ Model saved: {args.output_model}")
    model.wv.save_word2vec_format(args.output_bin, binary=True)
    print(f"✅ Binary vectors saved: {args.output_bin}")


#if __name__ == "__main__":
#    args = parse_args()
#    main(args)

##### USAGE FROM CONSOLE
from types import SimpleNamespace

args = SimpleNamespace(
    input_corpus='./data/pubmed_sample',
    input_mesh='./data/MeSH_graph.edgelist',
    input_dic='./data/MeSH_dic.pkl.gz',
    output_model='./pubmed_mesh_test',
    output_bin='./pubmed_mesh_test.bin',
    dimensions=200,
    walk_length=50,
    num_walks=2,
    window_size=5,
    epochs=5,
    min_count=5,
    sg=1,
    workers=8,
    p=2.0,
    q=1.0,
    directed=False
)

main(args)