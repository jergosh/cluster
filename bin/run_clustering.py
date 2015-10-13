import re
import sys
import glob
from os import path
from argparse import ArgumentParser
from subprocess import Popen

clustering_cmd = "bsub python bin/pdb_clustering.py --ensid {} --pdbmap data/pdb_map_Eutheria.tab --pdbfile data/pdb/pdb{}.ent --alnfile data/pdb_map/{}/{}.fa --mapfile data/pdb_map/{}/{}.tab --colorfile data/colorfiles/{}.tab --niter 10000 > {}"

argparser = ArgumentParser()
argparser.add_argument("--pdbmap", metavar="pdb_map", type=str, required=True)
argparser.add_argument("--outdir", metavar="dir", type=str, required=True)

args = argparser.parse_args()

for l in open(args.pdbmap):
    f = l.rstrip().split()
    ens, pdb = f[:2]
    ens_dir = ens[-2:]
    outfile = path.join(outdir, ens+'_'+pdb+'.dat')
    clustering = clustering_cmd.format(ens, pdb, ens_dir, ens, ens_dir, ens, ens, outfile)
    print >>sys.stderr, clustering

    p = Popen(clustering.split())
    p.wait()
