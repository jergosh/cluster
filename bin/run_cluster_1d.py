import os
from os import path
import pandas
import sys
import argparse
from subprocess import Popen
import glob 

import utils 

cluster1d_cmd = "python /hps/nobackup/goldman/gregs/cluster/bin/cluster_1d.py --infile {} --outfile {} --niter {}"

def main():
    argparser = argparse.ArgumentParser()

    argparser.add_argument('--indir', metavar='input_dir', type=str, required=True)
    argparser.add_argument('--filename', metavar='slr_file', type=str, default="ENS*.res")
    argparser.add_argument('--outdir', metavar='output_dir', type=str, required=True)
    argparser.add_argument('--logdir', metavar='log_dir', type=str, required=True)
    argparser.add_argument('--niter', metavar='niter', type=int, default=999)
    argparser.add_argument('--rerun', dest='rerun', action='store_true')
    argparser.add_argument('--queue', dest='queue', default="research", choices=["short", "research", "long"])

    argparser.set_defaults(rerun=False)


    args = argparser.parse_args()

    for f in glob.glob(path.join(args.indir, '*', args.filename)):
        basename = path.basename(f).rpartition('.')[0]
        basename_split = basename.split('_')
        out_id = '_'.join(basename_split[:3])
        subdir = basename_split[1][:2]
        outdir = path.join(args.outdir, subdir)
        logdir = path.join(args.logdir, subdir)
        outfile = path.abspath(path.join(outdir, out_id+'.res'))
        logfile = path.abspath(path.join(logdir, basename+'.log'))
        print outdir, out_id

        if args.rerun and  path.exists(outfile):
            print "Skipping", outfile
            continue

        utils.check_dir(outdir)
        utils.check_dir(logdir)
        cluster1d = cluster1d_cmd.format(path.abspath(f), outfile, args.niter)

        p = Popen([ 'bsub', '-R"affinity[core(1,same=socket,exclusive=(core, alljobs)): cpubind=core]"', '-q'+args.queue,
                    '-o'+ logfile, cluster1d ])
        p.wait()

if __name__ == "__main__":
    main()
