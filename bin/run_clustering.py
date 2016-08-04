import re
import sys
import glob
import csv
from os import path
from argparse import ArgumentParser
from subprocess import Popen

import numpy as np
import pandas

clustering_cmd = "python bin/pdb_clustering.py --pdbmap {} --pdbfile {} --outfile {} \
--thr {} --rerun_thr {} --stat {} {} --niter {} --rerun_iter {}"

argparser = ArgumentParser()
argparser.add_argument("--pdbmap", metavar="pdb_map", type=str, required=True)
argparser.add_argument("--pdbdir", metavar="pdb_dir", type=str, required=True)
argparser.add_argument("--outdir", metavar="out_dir", type=str, required=True)
argparser.add_argument("--logdir", metavar="log_dir", type=str, required=True)

argparser.add_argument("--method", metavar="method", type=str, choices=["cucala", "clumps"], required=True)
argparser.add_argument("--thr", metavar="thr", type=float, default=1.0)
argparser.add_argument("--rerun_thr", metavar="rerun_thr", type=float, default=0.001)
argparser.add_argument("--stat", metavar="thr", type=str, default="omega")
argparser.add_argument('--greater', dest='greater', action='store_true')
argparser.add_argument('--lesser', dest='greater', action='store_false')
argparser.set_defaults(greater=True)

argparser.add_argument("--niter", metavar="n_iter", type=int, required=False, default=0)
argparser.add_argument("--rerun_iter", metavar="rerun_iter", type=int, required=False, default=0)

args = argparser.parse_args()

def submit_clustering(df, pdbdir, thr, stat, greater, niter, rerun_thr, rerun_iter, outdir, logdir, method):
    # Write out a file -- stable_id - pdb_id - pdb_chain
    stable_id = df.stable_id.iloc[0]
    pdb_id = df.pdb_id.iloc[0]
    pdb_chain = df.pdb_chain.iloc[0]

    df_file = path.join(outdir, '_'.join([stable_id, pdb_id, pdb_chain])+'.tab')
    out_file = path.join(outdir, '_'.join([stable_id, pdb_id, pdb_chain])+'.out')
    log_file = path.join(logdir, '_'.join([stable_id, pdb_id, pdb_chain, str(thr)])+'.log')
    pdb_file = path.join(pdbdir, 'pdb'+pdb_id+'.ent')
    greater_val = "--greater" if greater else "--lesser"

    df.to_csv(df_file, sep="\t", quoting=csv.QUOTE_NONE)
    clustering = clustering_cmd.format(df_file,
                                       pdb_file,
                                       out_file,
                                       thr,
                                       rerun_thr,
                                       stat,
                                       greater_val,
                                       niter,
                                       rerun_iter,
                                       method)
    p = Popen([ 'bsub', '-o'+log_file, '-R"affinity[core(1,same=socket,exclusive=(core, alljobs)):cpubind=core]"', clustering ])
    p.wait()

    return df

pdb_map = pandas.read_table(args.pdbmap, dtype={ "stable_id": str, "pdb_id": str, "pdb_pos": str, "omega": np.float64 })
pdb_map.groupby(["stable_id", "pdb_id", "pdb_chain"]).apply(submit_clustering, args.pdbdir, args.thr, args.stat, args.greater, args.niter, args.rerun_thr, args.rerun_iter, args.outdir, args.logdir, args.method)
