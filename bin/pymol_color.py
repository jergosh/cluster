import os
from os import path
import glob
import sys
from operator import itemgetter
from argparse import ArgumentParser
import pandas as pd
# Include paths to 
# Is this the best way?
sys.path.extend(["/Library/Python/2.5/site-packages", 
                 "/Library/Python/2.5/site-packages/biopython-1.62-py2.5-macosx-10.9-x86_64.egg"])

from pymol import cmd
from Bio import SeqUtils
from Bio import AlignIO

argparser = ArgumentParser()

argparser.add_argument("--pdbstr", metavar="pdb_structure", type=str, required=True)
argparser.add_argument("--colorfile", metavar="color_file", type=str, required=False, default=None)

args = argparser.parse_args(os.environ['CMD'].split())

data_table = pd.read_table(open(args.colorfile), sep="\t")
# assert len(set(data_table.pdb_chain)) == 1
assert len(set(data_table.pdb_id)) == 1
pdb_name = data_table.pdb_id[0]
# cmd.load(path.join("data/pdb", 'pdb'+pdb_name+'.ent'))
cmd.load(args.pdbstr)

# This is weird.
seq_data = set()
# myspace = { 'f': lambda resi, resn: seq_data.add((int(resi), resn)) }
# cmd.iterate('(chain'+chain+')', "f(resi, resn)", space = myspace)

# resis = sorted(seq_data, key=itemgetter(0))
# seq = SeqUtils.seq1(''.join((c[1] for c in sorted(seq_data, key=itemgetter(0)))))

# for i, l in data_table.iterrows():
#     residue = str(l.pdb_pos)
#     if l.omega > 1.0:
#         r, g, b  = 1.0, 0.0, 0.0
#     else:
#         r, g, b  = 1.0, 1.0, 1.0

#     obj = '/pdb'+pdb_name+'//'+chain+'/'+residue
#     color = [r, g, b]

#     cmd.set_color("col"+residue, color)
#     cmd.color("col"+residue, obj)

def colour_pdb(df):
    print df

    chain = df.pdb_chain.iloc[0]

    # TODO: Where to get the chain from? Parse out of the table and assert that there is a 
    # single chain?
    cmd.select("chain"+chain, "chain "+chain)
    cmd.color("gray", "/"+pdb_name+'//'+chain)

    for i, l in data_table.iterrows():
        residue = str(l.pdb_pos)
        if l.omega == "grey":
            r, g, b  = 1.0, 1.0, 1.0
        elif l.omega == "red":
            r, g, b  = 1.0, 0.0, 0.0
        elif l.omega == "blue":
            r, g, b  = 0.0, 0.0, 1.0

        # obj = '/pdb'+pdb_name+'//'+chain+'/'+residue
        obj = '/'+pdb_name+'//'+chain+'/'+residue
        color = [r, g, b]

        cmd.set_color("col"+residue, color)
        cmd.color("col"+residue, obj)
    
    cmd.refresh()

data_table.groupby("pdb_chain").apply(colour_pdb)
# proj_dir = "/Users/greg/Documents/projects/slr_pipeline"
# ens, pdb_name, chain, score, length = None, None, None, None, None
# else:        
#     offset = 0
#     strsites_fn = path.join(proj_dir, 'data/ens/73/pdb_map/Eutheria', ens+'.tab')
#     strsites = open(strsites_fn)

#     strsites.readline() # Header
#     for l in strsites:
#         f = l.rstrip().split('\t')

#         if len(f) < 6:
#             continue
#         site, omega = int(f[2]), float(f[5])
#         obj = '/pdb'+pdb_name+'//'+chain+'/'+str(site+offset)

#         if omega > 1:
#             cmd.set_color("col"+str(site), [1, 0, 0])
#         else:
#             cmd.set_color("col"+str(site), [1, (1-omega)/1, (1-omega)/1])

#             cmd.color("col"+str(site), obj)
    
#             cmd.refresh()

