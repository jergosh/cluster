import os
from os import path
import pandas
import sys
import argparse
from subprocess import Popen
import glob 

cluster1d_cmd = "python /nfs/research2/goldman/gregs/cluster/bin/cluster_1d.py --infile {} --niter {}"

def main():
    argparser = argparse.ArgumentParser()

    argparser.add_argument('--indir', metavar='input_dir', type=str, required=True)
    argparser.add_argument('--filename', metavar='slr_file', type=str, default="ENS*.res")
    argparser.add_argument('--logdir', metavar='log_dir', type=str, required=True)
    argparser.add_argument('--niter', metavar='niter', type=int, default=999)

    args = argparser.parse_args()

    for f in glob.glob(path.join(args.indir, '*', args.filename)):
        basename = path.basename(f).rpartition('.')[0]
        print basename

        log_file = path.abspath(path.join(args.logdir, basename+'.log'))
        print path.abspath(f)
        cluster1d = cluster1d_cmd.format(path.abspath(f), args.niter)

        p = Popen([ 'bsub', '-o'+ log_file, cluster1d ])
        p.wait()

if __name__ == "__main__":
    main()
