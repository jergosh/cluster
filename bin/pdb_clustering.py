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
import cucala

import multiprocessing
import functools

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

def dist_min(r1, r2):
    dist = float('inf')

    for at1 in r1:
        for at2 in r2:
            d = at1 - at2

            if d < dist:
                dist = d

    return dist

# def signMWcont_iter(iter, coords, marks, dists):
#     marks_p = random.sample(marks, len(marks))
#     I, c, R, v =  cucala.MWcont(coords, marks_p, dists)

#     return I

# def signMWcont_multi(coords, marks, dists, niter, nthreads, pool):
#     maxI, maxCoords, maxR, maxV = cucala.MWcont(coords, marks, dists)

#     iter_partial = functools.partial(signMWcont_iter,
#                                      coords=coords, marks=marks, dists=dists)
    
#     results = pool.map(iter_partial, range(niter))
    
#     pval = float(1 + sum([I >= maxI for I in results ])) / (niter+1)
    
#     return maxI, maxCoords, maxR, maxV, pval

def cucala_pdb(sel_residues, all_residues, ids, dists, niter, pool):
    centroids = []
    marks = []

    for r in all_residues:
        # FIXME We're calculating the centroids twice
        centroids.append(centroid(r))
        
        marks.append(r in sel_residues)

    return cucala.signMWcont_multi(centroids, marks, ids, dists, niter, pool)

def run_cucala(sel_residues, all_residues, thr, niter, rerun_thr, rerun_iter, nthreads):
    rets = []
    centroids = [ centroid(r) for r in all_residues ]
    ids = [ r.id[0] + str(r.id[1]) + r.id[2] for r in all_residues ]

    # dists = cucala.order_dists(centroids)
    dists = cucala.order_dists(all_residues, dist_fun=dist_min)
    cluster_id = 1

    p = multiprocessing.Pool(nthreads)
    ret = cucala_pdb(sel_residues, all_residues, ids, dists, niter, p)
    # Output pre- and post-threshold p-values to separate files?
    if ret[4] < rerun_thr:
        print >>sys.stderr, ret[4], "rerunning..."
        ret = cucala_pdb(sel_residues, all_residues, ids, dists, rerun_iter, p)

    rets.append(ret)

    while ret[4] < thr:
        cluster_id += 1
        to_keep = [ i for i, item in enumerate(ids) if item not in ret[3] ]

        ids[:] = [ item for i, item in enumerate(ids) if i in to_keep ]
        all_residues[:] = [ item for i, item in enumerate(all_residues) if i in to_keep ]
        # = [ item for i, item in enumerate(all_residues) if i not in ret[1] ]
        # centroids = [ centroid(r) for r in all_residues ]
        centroids[:] = [ item for i, item in enumerate(centroids) if i in to_keep ]
        # dists = cucala.order_dists(centroids)
        dists = cucala.order_dists(all_residues, dist_fun=dist_min)
        
        ret = cucala_pdb(sel_residues, all_residues, ids, dists, niter, p)

        if ret[4] < rerun_thr:
            ret = cucala_pdb(sel_residues, all_residues, ids, dists, rerun_iter, p)

        if ret[4] < thr:
            rets.append(ret)

    return rets

# TODO Should have run_clumps() implemented as well
def run_clumps(sel_residues, all_residues, thr, niter, rerun_thr, rerun_iter):
    WAP = clumps(sel_residues)

    pval = run_sim(pdb_chain, len(sel_residues), WAP, niter)
    print >>sys.stderr, pval,
    if pval < rerun_thr:
        pval = run_sim(pdb_chain, len(sel_residues), WAP, rerun_iter)
        print >>sys.stderr, pval
    else:
        print >>sys.stderr

    return pval
        
def clumps(residues, thr=6):
    WAP = 0.0
    centroids = [ centroid(r) for r in residues ]

    for i, r1 in enumerate(centroids):
        for r2 in centroids[i+1:]:
            d = dist(r1, r2)
            WAP += exp(- d**2 / (2*thr**2))

    return WAP

def run_sim(pdb_chain, n_residues, score, niter, method):
    sim_WAP = []
    for i in range(niter):
        sim_residues = random.sample(pdb_chain, n_residues)
        sim_WAP.append(method(sim_residues))

    pval = sum([ score < w for w in sim_WAP ]) * 1.0 / (niter+1)
    return pval

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

def process_pdb(df, pdbfile, thr, stat, greater, niter, rerun_thr, rerun_iter, outfile, method, sign_thr, nthreads):
    pdb_id = df.pdb_id.iloc[0]
    stable_id = df.stable_id.iloc[0]
    chain_id = df.pdb_chain.iloc[0]

    if greater:
        op = operator.gt
    else:
        op = operator.lt

    try:
        pdb = p.get_structure(pdb_id, pdbfile)
        pdb_chain = pdb[0][chain_id]
    except IOError, e:
        print >>sys.stderr, "PDB file", pdb_id, "missing!"
        return df

    all_residues = []
    sel_residues = []
    sel_coords = []
    for i, row in df.iterrows():
        res_id = parse_coord(row.pdb_pos)
        try:
            r = pdb_chain[res_id]
        except KeyError, e:
            r = find_sequential(pdb_chain, res_id)
            if r is None:
                raise e

            print >>sys.stderr, r.get_id()

        all_residues.append(r)

        if op(row[stat], thr):
            sel_coords.append(r['CA'].get_coord())
            sel_residues.append(r)

    print >>sys.stderr, stable_id, pdb_id, chain_id, '('+str(len(sel_coords))+')',

    if len(sel_residues) < 2:
        return df
    if method == "clumps":
        pval = run_clumps(sel_residues, all_residues, sign_thr, niter, rerun_thr, rerun_iter)
        print >>outfile, '\t'.join([ str(it) for it in 
                                     [ stable_id, pdb_id, pdb_chain.id, len(pdb_chain), len(all_residues), '[]', pval ] ])

    elif method == "cucala":
        rets = run_cucala(sel_residues, all_residues, sign_thr, niter, rerun_thr, rerun_iter, nthreads)

        for ret in rets:
            print >>outfile, '\t'.join([ str(i) for i in [ stable_id, pdb_id, pdb_chain.id, len(pdb_chain), len(all_residues), ret[3], ret[4] ] ])

    # print '\t'.join([ str(it) for it in 
    #                   [ cath_id, pdb_id, len(pdb_chain), len(residues), pval ] ])
    outfile.close()
    return df

p = PDB.PDBParser(QUIET=True)

argparser = ArgumentParser()
argparser.add_argument("--pdbmap", metavar="pdb_map", type=str, required=True)
argparser.add_argument("--pdbfile", metavar="pdb_dir", type=str, required=True)
argparser.add_argument('--outfile', metavar='out_file', type=str, required=True)

argparser.add_argument("--method", metavar="method", type=str, choices=["cucala", "clumps"], required=True)
argparser.add_argument("--sign_thr", metavar="sign_thr", type=float, default=0.05)
argparser.add_argument("--thr", metavar="thr", type=float, default=1.0)
argparser.add_argument("--rerun_thr", metavar="rerun_thr", type=float, default=0.001)
argparser.add_argument("--stat", metavar="thr", type=str, default="omega")
argparser.add_argument('--greater', dest='greater', action='store_true')
argparser.add_argument('--lesser', dest='greater', action='store_false')
argparser.set_defaults(greater=True)

argparser.add_argument("--niter", metavar="n_iter", type=int, required=False, default=0)
argparser.add_argument("--rerun_iter", metavar="rerun_iter", type=int, required=False, default=0)
# argparser.add_argument("--preference", metavar="preference", type=int, required=False, default=-50)

argparser.add_argument("--nthreads", metavar="n_threads", type=int, required=False, default=1)

args = argparser.parse_args()

pdb_map = pandas.read_table(args.pdbmap, dtype={ "stable_id": str, "pdb_id": str, "pdb_pos": str, "omega": np.float64 })

# pdb_map = pandas.read_table(args.pdbmap, 
#                             names=[ "ens_id", "ens_coord", "uniprot_id", "uniprot_coord", 
#                                       "pdb_id", "pdb_chain", "pdb_coord", "SS", "rsa", "omega" ],
#                             dtype={ "ens_id": str, "ens_coord": int, "uniprot_id": str, "uniprot_coord": int,
#                                     "pdb_id": str, "pdb_chain": str, "pdb_coord": str, "SS": str,
#                                     "rsa": np.float64, "omega": np.float64 })

# pdb_map.groupby(["cath_id", "pdb_id"]).apply(process_pdb, args.pdbfile)
outfile = open(args.outfile, 'w')
# pdb_map.groupby(["stable_id", "pdb_id", "pdb_chain"]).apply(process_pdb, args.pdbdir, args.thr, args.stat, args.greater, args.niter, args.rerun_thr, args.rerun_iter, outfile)
process_pdb(pdb_map, args.pdbfile, args.thr, args.stat, args.greater, args.niter, args.rerun_thr, args.rerun_iter, outfile, args.method, args.sign_thr, args.nthreads)
