import os
from os import path
import glob
import sys
from operator import itemgetter
from argparse import ArgumentParser
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


proj_dir = "/Users/greg/Documents/projects/slr_pipeline"
pdb_map_file = open(path.join(proj_dir, "data/pdb_map_Eutheria.tab"))

ens, pdb_name, chain, score, length = None, None, None, None, None
for l in pdb_map_file:
    ens, pdb_name, chain, score, length = l.rstrip().split('\t')[:5]
    if pdb_name == args.pdbstr:
        break
else:
    raise Exception

cmd.load(path.join(proj_dir, "data/pdb", 'pdb'+pdb_name+'.ent'))
cmd.select("chain"+chain, "chain "+chain)

# This is weird.
seq_data = set()
myspace = { 'f': lambda resi, resn: seq_data.add((int(resi), resn)) }
cmd.iterate('(chain'+chain+')', "f(resi, resn)", space = myspace)

resis = sorted(seq_data, key=itemgetter(0))
seq = SeqUtils.seq1(''.join((c[1] for c in sorted(seq_data, key=itemgetter(0)))))

offset = 0
cmd.color("bluewhite", "/pdb"+pdb_name+'//'+chain)

if args.colorfile is not None:
    colorfile = open(args.colorfile)
    
    for l in colorfile:
        f = l.rstrip().split('\t')
        residue = f[0].strip()
        r, g, b  = (float(i) for i in f[1:])
        obj = '/pdb'+pdb_name+'//'+chain+'/'+residue
        color = [r/255, g/255, b/255]

        cmd.set_color("col"+residue, color)
        cmd.color("col"+residue, obj)
    
    cmd.refresh()
else:        
    strsites_fn = path.join(proj_dir, 'data/ens/73/pdb_map/Eutheria', ens+'.tab')
    strsites = open(strsites_fn)

    strsites.readline() # Header
    for l in strsites:
        f = l.rstrip().split('\t')

        if len(f) < 6:
            continue
        site, omega = int(f[2]), float(f[5])
        obj = '/pdb'+pdb_name+'//'+chain+'/'+str(site+offset)

        if omega > 1:
            cmd.set_color("col"+str(site), [1, 0, 0])
        else:
            cmd.set_color("col"+str(site), [1, (1-omega)/1, (1-omega)/1])

            cmd.color("col"+str(site), obj)
    
            cmd.refresh()
