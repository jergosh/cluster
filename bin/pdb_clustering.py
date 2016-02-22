import re
import sys
import glob
from os import path
from argparse import ArgumentParser
import random
# import utils
from collections import defaultdict
from pprint import pprint
from array import array
from math import sqrt, exp
import operator

from Bio import SeqIO
from Bio import AlignIO
from Bio import SeqUtils
from Bio import PDB
import pandas
import numpy as np

import networkx as nx

re_resid = re.compile("(-?[0-9]+)([A-Z]*)")

def isnan(f):
    return not (f == f)

def parse_coord(coord):
    n, ic = re_resid.match(coord).groups()
    if ic == '':
        ic = ' '

    return ' ', int(n), ic

# TODO Put back the slow version and compare results?

def compute_neighbours(chain, dist_thr):
    neighbour_map = defaultdict(set)

    for r1 in chain:
        if SeqUtils.seq1(r1.resname) == 'X':
            continue
        
        for at1 in r1:
            for r2 in chain:
                if SeqUtils.seq1(r2.resname) == 'X':
                    continue
                    
                for at2 in r2:
                    if at1 - at2 < dist_thr:
                        neighbour_map[r1].add(r2)
                        neighbour_map[r2].add(r1)
                        continue

    return neighbour_map
        
def rcn_clustering(neighbour_map, chain, sign_sites):
    G = nx.Graph()
    for r in sign_sites:
        G.add_node(r.id)
        for nr in neighbour_map[r]:
            if nr in sign_sites:
                G.add_edge(r.id, nr.id)

    ccs = nx.algorithms.components.connected_components(G)
    
    res_map = {}
    [ [ res_map.__setitem__(r, i) for r in l ] for i, l in enumerate(ccs) ]

    labels = [ res_map[r.id] for r in sign_sites ]
    return labels

def centroid(residue):
    c = array('f', [0, 0, 0])

    k = 0
    for at in residue:
        c += at.get_coord()
        k += 1
        
    return c / k 

def dist(at1, at2):
    return sqrt(reduce(operator.add, [ (c[0]-c[1])**2 for c in zip(at1, at2) ]))
        
def clumps(residues, thr=6):
    WAP = 0.0
    centroids = [ centroid(r) for r in residues ]

    for i, r1 in enumerate(centroids):
        for r2 in centroids[i+1:]:
            d = dist(r1, r2)
            WAP += exp(- d**2 / (2*thr**2))

    return WAP

# TODO: Is it worth it to make this into a class (or refactor altogether)
d = lambda: defaultdict(int)
class ClusterList():
    def __init__(self):
        self.cl_found = defaultdict(d)
        self.cl_max = defaultdict(int)
        self.cl_number = defaultdict(int)

    def collapse_clusters(self, labels):
        clusters = defaultdict(int)
        cluster_sizes = defaultdict(int)
        for l in labels:
            if isnan(l):
                print >>sys.stderr, labels
                return {}
            clusters[l] += 1

        for cl_size in clusters.values():
            cluster_sizes[cl_size] += 1

        return cluster_sizes

    def process_clusters(self, labels):
        cluster_sizes = self.collapse_clusters(labels)
        self.cl_max[max(cluster_sizes.keys())] += 1
        self.cl_number[max(labels)+1] += 1

        for cl_size in cluster_sizes:
            self.cl_found[cl_size][cluster_sizes[cl_size]] += 1

    def get_pvals(self, niter):
        for n in sorted(self.cl_found, reverse=True):
            for k in sorted(self.cl_found[n], reverse=True):
                if n > 1:
                    self.cl_found[n-1][k] += self.cl_found[n][k]
                if k > 1:
                    self.cl_found[n][k-1] += self.cl_found[n][k]

                self.cl_found[n][k] = min(self.cl_found[n][k]/float(niter), 1.0)

        # Biggest cluster size
        for k in sorted(self.cl_max, reverse=True):
            if k > 1:
                self.cl_max[k-1] += self.cl_max[k]

            self.cl_max[k] = min(self.cl_max[k]/float(niter), 1.0)

        # Number of clusters
        for k in range(1, max(self.cl_number)+1):
            if k < max(self.cl_number):
                self.cl_number[k+1] += self.cl_number[k]

            self.cl_number[k] = min(self.cl_number[k]/float(niter), 1.0)

def find_sequential(chain, res_id):
    for r in chain:
        if r.get_id()[1:] == res_id[1:]:
            return r

    return None

def process_pdb(df, pdbdir):
    cath_id = df.cath_id.iloc[0]
    pdb_id = df.pdb_id.iloc[0]
    # ens_id = df.ens_id.iloc[0]
    # chain_id = df.pdb_chain.iloc[0]
    # print >>sys.stderr, ens_id, pdb_id, chain_id
    print >>sys.stderr, pdb_id
    try:
        # pdb = p.get_structure(pdb_id, path.join(pdbdir, 'pdb'+pdb_id+'.ent'))
        pdb = p.get_structure(pdb_id, path.join(pdbdir, pdb_id))
        # pdb_chain = pdb[0][chain_id]
        pdb_chain = list(pdb[0])[0]
    except IOError, e:
        print >>sys.stderr, "PDB file", pdb_id, "missing!"
        return

    r_coords = []
    residues = []
    for i, row in df.iterrows():
        if row.omega > args.omega_thr:
            res_id = parse_coord(row.pdb_coord)
            try:
                r = pdb_chain[res_id]
            except KeyError, e:
                r = find_sequential(pdb_chain, res_id)
                if r is None:
                    raise e
                print >>sys.stderr, r.get_id()

            r_coords.append(r['CA'].get_coord())
            residues.append(r)

    if len(r_coords) < 2:
        # print >>sys.stderr, "<= 1 residue with omega > {}. Exiting".format(args.omega_thr)
        return

    WAP = clumps(residues)
    #neighbour_map = compute_neighbours(pdb_chain, args.dist_thr)
    #labels = rcn_clustering(neighbour_map, pdb_chain, residues)

    # print >>sys.stderr, "Running simulations...",
    sim_WAP = []
    for i in range(args.niter):
        sim_residues = random.sample(pdb_chain, len(residues))
        sim_WAP.append(clumps(sim_residues))

    # print WAP, sum(sim_WAP) / len(sim_WAP)
    pval = sum([ WAP < w for w in sim_WAP ]) * 1.0 / args.niter
    print '\t'.join([ str(it) for it in 
                      [ cath_id, pdb_id, len(pdb_chain), len(residues), pval ] ])
    # print '\t'.join([ str(it) for it in 
    #                   [ ens_id, pdb_id, len(pdb_chain), len(residues), pval ] ])

    # clusterlist = ClusterList()
    # i = 0
    # while i < args.niter:
    #     sim_residues = random.sample(pdb_chain, len(residues)) 
    #     sim_labels = rcn_clustering(neighbour_map, pdb_chain, sim_residues)

    #     if any([ isnan(it) for it in sim_labels ]):
    #         continue

    #     clusterlist.process_clusters(sim_labels)
    #     i += 1
    

    # print >>sys.stderr, "done."
    # clusterlist.get_pvals(args.niter)

    # clusters = clusterlist.collapse_clusters(labels)
    # combined_pval = 1.0
    # for n in sorted(clusters):
    #     k = clusters[n]
    #     combined_pval *= clusterlist.cl_found[n][k]

    # print '\t'.join([ str(it) for it in 
    #                   [ ens_id, pdb_id, len(pdb_chain), len(labels), combined_pval, clusterlist.cl_max[max(clusters.keys())], clusterlist.cl_number[max(labels)+1] ] ])

p = PDB.PDBParser(QUIET=True)

argparser = ArgumentParser()
argparser.add_argument("--pdbmap", metavar="pdb_map", type=str, required=True)
argparser.add_argument("--pdbdir", metavar="pdb_dir", type=str, required=True)
argparser.add_argument("--dist_thr", metavar="thr", type=float, default=4.0)
argparser.add_argument("--omega_thr", metavar="thr", type=float, default=1.0)

argparser.add_argument("--niter", metavar="n_iter", type=int, required=False, default=0)
argparser.add_argument("--preference", metavar="preference", type=int, required=False, default=-50)


args = argparser.parse_args()

pdb_map = pandas.read_table(args.pdbmap, dtype={ "cath_id": str, "pdb_id": str, "pdb_coord": str, "omega": np.float64 })

# pdb_map = pandas.read_table(args.pdbmap, 
#                             names=[ "ens_id", "ens_coord", "uniprot_id", "uniprot_coord", 
#                                       "pdb_id", "pdb_chain", "pdb_coord", "SS", "rsa", "omega" ],
#                             dtype={ "ens_id": str, "ens_coord": int, "uniprot_id": str, "uniprot_coord": int,
#                                     "pdb_id": str, "pdb_chain": str, "pdb_coord": str, "SS": str,
#                                     "rsa": np.float64, "omega": np.float64 })

pdb_map.groupby(["cath_id"]).apply(process_pdb, args.pdbdir)
# pdb_map.groupby(["ens_id", "pdb_id", "pdb_chain"]).apply(process_pdb, args.pdbdir)
# process_pdb(pdb_map)
