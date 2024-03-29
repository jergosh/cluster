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

def rid2str(r):
    return r.id[0] + str(r.id[1]) + r.id[2]

## My graph-based methods
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

def run_graph(sel_residues, all_residues, thr, niter, rerun_thr, rerun_iter, dist_thr=6):
    neighbour_map = compute_neighbours(all_residues, dist_thr=dist_thr)

    max_clust, max_clust_id, n_clusts, res_map = graph_clustering(neighbour_map, sel_residues)
    max_labels = [ r[0] + str(r[1]) + r[2] for r, i, in res_map.items() if i == max_clust_id ]
    # labels = [ res_map[r.id] for r in sel_residues ]

    pval_max, pval_n = graph_sim(sel_residues, all_residues, neighbour_map, max_clust, n_clusts, niter)
    
    if pval_max < rerun_thr or pval_n < rerun_thr:
        pval_max, pval_n = graph_sim(sel_residues, all_residues, neighbour_map, max_clust, n_clusts, rerun_iter)

    return [[pval_max, max_labels], [pval_n]]

def graph_sim(sel_residues, all_residues, neighbour_map, max_clust, n_clusts, niter):
    pval_max, pval_n = 1.0, 1.0

    for i in range(niter):
        sim_residues = random.sample(all_residues, len(sel_residues))
        max_clust_sim, max_clust_id_sim, n_clusts_sim, sim_map = graph_clustering(neighbour_map, sim_residues)

        if max_clust_sim >= max_clust:
            pval_max += 1

        if n_clusts_sim <= n_clusts:
            pval_n += 1

    pval_max = pval_max / (niter+1)
    pval_n = pval_n / (niter+1)

    return pval_max, pval_n


def graph_clustering(neighbour_map, sel_residues):
    G = nx.Graph()
    for r in sel_residues:
        G.add_node(r.id)
        for nr in neighbour_map[r]:
            if nr in sel_residues:
                G.add_edge(r.id, nr.id)

    ccs = list(nx.algorithms.components.connected_components(G))
    max_clust, max_clust_id = sorted([ (len(cc), i) for i, cc in enumerate(ccs) ], key=operator.itemgetter(0))[-1]
    n_clusts = len(ccs)
    
    res_map = {}
    [ [ res_map.__setitem__(r, i) for r in l ] for i, l in enumerate(ccs) ]

    return max_clust, max_clust_id, n_clusts, res_map

## Cucala
def centroid(residue):
    c = array('f', [0, 0, 0])

    k = 0
    for at in residue:
        c += at.get_coord()
        k += 1
        
    return c / k 

# def dist(at1, at2):
#     return sqrt(reduce(operator.add, [ (c[0]-c[1])**2 for c in zip(at1, at2) ]))

def dist(at1, at2):
    return sqrt((at1[0]-at2[0])**2 + (at1[1]-at2[1])**2 + (at1[2]-at2[2])**2)

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

    # return cucala.signMWcont_multi(centroids, marks, ids, dists, niter, pool)
    return cucala.signMWcont(centroids, marks, ids, dists, niter)

def run_cucala(sel_residues, all_residues, mode, thr, niter, rerun_thr, rerun_iter, nthreads):
    rets = []
    centroids = [ centroid(r) for r in all_residues ]
    ids = [ rid2str(r) for r in all_residues ]

    dists = cucala.order_dists(centroids)
    cluster_id = 1

    # p = multiprocessing.Pool(nthreads, maxtasksperchild=1000)
    p = None
    ret = cucala_pdb(sel_residues, all_residues, ids, dists, niter, p)
    # Output pre- and post-threshold p-values to separate files?
    
    if ret[4] < rerun_thr:
        print >>sys.stderr, ret[4], "rerunning..."
        ret = cucala_pdb(sel_residues, all_residues, ids, dists, rerun_iter, p)

    rets.append(ret)

    while ret[4] < thr:
        cluster_id += 1

        to_keep = [ i for i, item in enumerate(ids) if item not in ret[3] ]
        print >>sys.stderr, to_keep

        ids[:] = [ item for i, item in enumerate(ids) if i in to_keep ]
        all_residues[:] = [ item for i, item in enumerate(all_residues) if i in to_keep ]
        # = [ item for i, item in enumerate(all_residues) if i not in ret[1] ]
        # centroids = [ centroid(r) for r in all_residues ]
        centroids[:] = [ item for i, item in enumerate(centroids) if i in to_keep ]
        dists = cucala.order_dists(centroids)
        
        ret = cucala_pdb(sel_residues, all_residues, ids, dists, niter, p)
        
        if ret[4] < rerun_thr:
            print >>sys.stderr, ret[4], "rerunning..."
            ret = cucala_pdb(sel_residues, all_residues, ids, dists, rerun_iter, p)

        if ret[4] < thr:
            rets.append(ret)

    return rets


# CLUMPS
def run_clumps(sel_residues, all_residues, thr, niter, rerun_thr, rerun_iter):
    WAP = clumps(sel_residues)

    pval = run_sim(all_residues, len(sel_residues), WAP, niter)

    print >>sys.stderr, pval,
    if pval < rerun_thr:
        print >>sys.stderr, pval
    else:
        print >>sys.stderr
        
    print clumps(sel_residues)
    return pval
        
def clumps(residues, thr=6):
    WAP = 0.0
    # coords = [ centroid(r) for r in residues ]
    coords = [ r['CA'].get_coord() for r in residues ]

    for i, r1 in enumerate(residues):
        for r2 in residues[i+1:]:
            d = dist(r1['CA'].get_coord(), r2['CA'].get_coord())
            WAP += exp(- d**2 / (2*thr**2))

    return WAP

def run_sim(pdb_chain, n_residues, score, niter):
    sim_WAP = []
    
    for i in range(niter):
        sim_residues = random.sample(pdb_chain, n_residues)
        w = clumps(sim_residues)

        sim_WAP.append(w)

    pval = (1+sum([ score < w for w in sim_WAP ])) * 1.0 / (niter+1)
    return pval


def find_sequential(chain, res_id):
    for r in chain:
        if r.get_id()[1:] == res_id[1:]:
            return r

    return None

def process_pdb(df, pdbfile, thr, stat, greater, niter, rerun_thr, rerun_iter, dist_thr, outfile, method, mode, sign_thr, nthreads):
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

        if op(row[stat], thr) and row["omega"] > 1.0:
            sel_coords.append(r['CA'].get_coord())
            sel_residues.append(r)

    print >>sys.stderr, stable_id, pdb_id, chain_id, '('+str(len(sel_coords))+')',

    if len(sel_residues) < 2:
        return df
    if method == "clumps":
        pval = run_clumps(sel_residues, all_residues, sign_thr, niter, rerun_thr, rerun_iter)
        ids = [ rid2str(r) for r in sel_residues ]
        print >>outfile, '\t'.join([ str(it) for it in 
                                     [ stable_id, pdb_id, pdb_chain.id, len(pdb_chain), len(all_residues), ids, pval ] ])

    elif method == "cucala":
        rets = run_cucala(sel_residues, all_residues, mode, sign_thr, niter, rerun_thr, rerun_iter, nthreads)

        for ret in rets:
            print >>outfile, '\t'.join([ str(i) for i in [ stable_id, pdb_id, pdb_chain.id, len(pdb_chain), len(all_residues), ret[3], ret[4] ] ])

    # print '\t'.join([ str(it) for it in 
    #                   [ cath_id, pdb_id, len(pdb_chain), len(residues), pval ] ])
    elif method == "gr":
        rets = run_graph(sel_residues, all_residues, sign_thr, niter, rerun_thr, rerun_iter, dist_thr)

        ret_max = rets[0]
        ret_n = rets[1]
        print >>outfile, '\t'.join([ str(i) for i in [ stable_id, pdb_id, pdb_chain.id, len(pdb_chain), len(all_residues), "gr_max", rets[0][1], rets[0][0]  ] ])
        print >>outfile, '\t'.join([ str(i) for i in [ stable_id, pdb_id, pdb_chain.id, len(pdb_chain), len(all_residues), "gr_n", "[]", rets[1][0]  ] ])
    
    outfile.close()
    return df

p = PDB.PDBParser(QUIET=True)

argparser = ArgumentParser()
argparser.add_argument("--pdbmap", metavar="pdb_map", type=str, required=True)
argparser.add_argument("--pdbfile", metavar="pdb_dir", type=str, required=True)
argparser.add_argument('--outfile', metavar='out_file', type=str, required=True)

argparser.add_argument("--method", metavar="method", type=str, choices=["cucala", "clumps", "gr"], required=True)
argparser.add_argument("--mode", metavar="mode", type=str, choices=["discrete", "continuous"], default="discrete")
argparser.add_argument("--sign_thr", metavar="sign_thr", type=float, default=0.05)
argparser.add_argument("--thr", metavar="thr", type=float, default=0.05)
argparser.add_argument("--rerun_thr", metavar="rerun_thr", type=float, default=0.001)
argparser.add_argument("--dist_thr", metavar="dist_thr", type=float, default=6)
argparser.add_argument("--stat", metavar="stat", type=str, default="Adj.Pval")
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
process_pdb(pdb_map, args.pdbfile, args.thr, args.stat, args.greater, args.niter, args.rerun_thr, args.dist_thr, args.rerun_iter, outfile, args.method, args.mode, args.sign_thr, args.nthreads)
