import os
from os import path
import pandas
import sys
import argparse
from subprocess import Popen
import glob 

cluster1d_cmd = "python /hps/nobackup/goldman/gregs/cluster/bin/cluster_1d.py --infile {} --outfile {} --niter {}"

def main():
    argparser = argparse.ArgumentParser()

    argparser.add_argument('--indir', metavar='input_dir', type=str, required=True)
    argparser.add_argument('--filename', metavar='slr_file', type=str, default="ENS*.res")
    argparser.add_argument('--outdir', metavar='output_dir', type=str, required=True)
    argparser.add_argument('--logdir', metavar='log_dir', type=str, required=True)
    argparser.add_argument('--niter', metavar='niter', type=int, default=999)

    args = argparser.parse_args()

    for f in glob.glob(path.join(args.indir, '*', args.filename)):
        basename = path.basename(f).rpartition('.')[0]
        out_id = '_'.join(basename.split('_')[:3])
        print out_id
        outfile = path.abspath(path.join(args.outdir, 'out_id'+'.res'))

        logfile = path.abspath(path.join(args.logdir, basename+'.log'))
        cluster1d = cluster1d_cmd.format(path.abspath(f), outfile, args.niter)

        p = Popen([ 'bsub', '-o'+ logfile, cluster1d ])
        p.wait()

if __name__ == "__main__":
    main()
