import os
from os import path
import sys
import argparse
import glob

def main():
    argparser = argparse.ArgumentParser()

    argparser.add_argument('--indir', metavar='slr_file', type=str, required=True)
    argparser.add_argument('--outfile', metavar='out_file', type=str, required=True)

    args = argparser.parse_args()

    outfile = open(args.outfile, 'w')
    
    for f in glob.glob(path.join(args.indir, '*', '*.res')):
        basename = path.basename(f)

        ens_id = basename.split('_')[0]
        dataset = '_'.join(basename.split('_')[1:3])
        print ens_id

        for l in open(f):
            ## TERRIBLE
            f = eval(l)
            print >>outfile, '\t'.join([ ens_id, dataset ] + f)


if __name__ == "__main__":
    main()

