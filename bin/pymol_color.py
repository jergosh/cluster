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
argparser.add_argument("--session", metavar="session_name", type=str, default=None)

args = argparser.parse_args(os.environ['CMD'].split())

# Could be more explicit about types here
data_table = pd.read_table(open(args.colorfile), sep="\t", dtype={ "pdb_pos": str })
print data_table
assert len(set(data_table.pdb_id)) == 1
pdb_name = data_table.pdb_id.iloc[0]
cmd.load(args.pdbstr)

# This assumes there is a single PDB chain (currently the input df is split on pdb chain)
def colour_pdb(df):
    chain = str(df.pdb_chain.iloc[0])

    # cmd.select("chain"+chain, "chain "+chain)
    cmd.hide("lines", "chain "+chain)
    cmd.show("cartoon", "chain "+chain)
    cmd.color("gray90", "/pdb"+pdb_name+'//'+chain)

    for i, l in data_table.iterrows():
        residue = str(l.pdb_pos)
        obj = '/pdb'+pdb_name+'//'+chain+'/'+residue

        cmd.color(l.color, obj)

def make_sels(df):
    group = df.group.iloc[0]
    chain = df.pdb_chain.iloc[0]
    
    resnames = '+'.join(df.pdb_pos)
    print resnames

    cmd.select(group, "chain "+chain+" and resi "+resnames )


# cmd.hide("everything")
cmd.show("sticks", "hetatm")
data_table.groupby("pdb_chain").apply(colour_pdb)
data_table.groupby("group").apply(make_sels)
cmd.refresh()
    
if session is not None:
    print "Saving to", args.session
    cmd.save(args.session)

