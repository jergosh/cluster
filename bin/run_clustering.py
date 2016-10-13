import re
import sys
import glob
import csv
from os import path
from argparse import ArgumentParser
from subprocess import Popen
import operator

import numpy as np
import pandas

clustering_cmd = "python bin/pdb_clustering.py --pdbmap {} --pdbfile {} --outfile {} \
--thr {} --rerun_thr {} --stat {} {} --niter {} --rerun_iter {} --method {} --sign_thr {} --nthreads {}"

argparser = ArgumentParser()
argparser.add_argument("--pdbmap", metavar="pdb_map", type=str, required=True)
argparser.add_argument("--pdbdir", metavar="pdb_dir", type=str, required=True)
argparser.add_argument("--outdir", metavar="out_dir", type=str, required=True)
argparser.add_argument("--logdir", metavar="log_dir", type=str, required=True)

argparser.add_argument("--method", metavar="method", type=str, choices=["cucala", "clumps", "gr"], required=True)
argparser.add_argument("--sign_thr", metavar="sign_thr", type=float, default=0.05)
argparser.add_argument("--rerun_thr", metavar="rerun_thr", type=float, default=0.001)
argparser.add_argument("--stat", metavar="stat", type=str, default="Adj.Pval")
argparser.add_argument("--thr", metavar="thr", type=float, default=0.05)
argparser.add_argument('--greater', dest='greater', action='store_true')
argparser.add_argument('--lesser', dest='greater', action='store_false')
argparser.set_defaults(greater=True)

argparser.add_argument("--niter", metavar="n_iter", type=int, required=False, default=0)
argparser.add_argument("--rerun_iter", metavar="rerun_iter", type=int, required=False, default=0)

argparser.add_argument("--nthreads", metavar="n_threads", type=int, required=False, default=1)

args = argparser.parse_args()

def submit_clustering(df, pdbdir, thr, stat, greater, niter, rerun_thr, rerun_iter, outdir, logdir, method, sign_thr, nthreads):
    stable_id = df.stable_id.iloc[0]
    pdb_id = df.pdb_id.iloc[0]
    pdb_chain = df.pdb_chain.iloc[0]

    # Check if there are enough sites
    if greater:
        op = operator.gt
    else:
        op = operator.lt

    n_check = 0
    for i, row in df.iterrows():
        if op(row[stat], thr) and row["omega"] > 1.0:
            n_check += 1

    if n_check < 2:
        print >>sys.stderr, "Skipping", stable_id, pdb_id
        return df

    # Write out a file -- stable_id - pdb_id - pdb_chain
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
                                       method,
                                       sign_thr,
                                       nthreads)
    # p = Popen([ 'bsub', '-o'+log_file, '-n'+str(nthreads), '-R"affinity[core({},same=socket,exclusive=(core, alljobs)):cpubind=core]"'.format(nthreads), clustering ])
    p = Popen([ 'bsub', '-M8192', '-o'+log_file, '-n'+str(nthreads), clustering ])
    p.wait()

    return df

pdb_map = pandas.read_table(args.pdbmap, dtype={ "stable_id": str, "pdb_id": str, "pdb_pos": str, "omega": np.float64 })
pdb_map.groupby(["stable_id", "pdb_id", "pdb_chain"]).apply(submit_clustering, args.pdbdir, args.thr, args.stat, args.greater, args.niter, args.rerun_thr, args.rerun_iter, args.outdir, args.logdir, args.method, args.sign_thr, args.nthreads)
