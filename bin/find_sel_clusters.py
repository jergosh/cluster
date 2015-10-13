from os import path
import glob
import warnings
# Suppress warnings from PDBParser
warnings.simplefilter("ignore")
from multiprocessing import Pool

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio import SeqUtils
from Bio import PDB
from Bio import AlignIO
from Bio.Align.Applications import TCoffeeCommandline

import networkx as nx

# Distance threshold
DIST_THR = 8

aln_dir = "data/mammals_slr_dataset/web"
pdb_dir = "data/pdb"
sites_dir = "data/strsites"
p = PDB.PDBParser()

best_pdb_file = open("results/best_pdb_biomart.txt")

def process_line(l):
    ens, pdb_name, chain, score, length = l.rstrip().split('\t')[0:5]
    # print pdb_name, ens
    pdb_fn = path.join(pdb_dir, 'pdb' + pdb_name + '.ent')
    pdb = p.get_structure(ens, pdb_fn)

    # The sequence from the PDB file
    chain = pdb[0][chain]
    # chain_seq = Seq(SeqUtils.seq1(''.join([ r.resname for r in chain ])))
    # Remove the Xs
    # chain_seq = chain_seq.ungap('X')

    # FIXME: Do we still need this offset?
    # offset = iter(chain).next().id[1]

    sites_fn = path.join(sites_dir, ens+'_strsites.csv')
    if not path.exists(sites_fn):
        return
    sites = open(sites_fn)
    nsites = 0
    sign_sites = set()
    for l in sites:
        site, omega = l.rstrip().split('\t')
        site, omega = int(site), float(omega)
        if omega > 1:
            # FIXME sign_sites.add(site+offset)
            sign_sites.add(site)
            nsites += 1

    if nsites < 2:
        return

    # The 'residue contact network'
    G = nx.Graph()
    for r1 in chain:
        if SeqUtils.seq1(r1.resname) == 'X':
            continue

        for r2 in chain:
            if SeqUtils.seq1(r2.resname) == 'X':
                continue

            if r1 == r2:
                continue
            if r1.id[1] in sign_sites and r2.id[1] in sign_sites:
                ca1 = r1['CA']
                ca2 = r2['CA']

                if ca1 - ca2 < DIST_THR:
                    G.add_edge(r1.id, r2.id)

    ccs = nx.algorithms.components.connected_components(G)
    total_size = len([ item for sublist in ccs if len(sublist) > 3 for item in sublist ])

    if total_size:
        print pdb_name, ens, total_size, len([ item for sublist in ccs for item in sublist ])

pool = Pool(7)
pool.map(process_line, best_pdb_file.readlines())
# pool.map(process_line, ["ENSP00000312741	2zv2	A	Q96RR4		8	298	159	449	158	448"])

