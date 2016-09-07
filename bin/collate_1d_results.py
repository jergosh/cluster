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
        dataset = '_'.join(basename.partition('.')[0].split('_')[1:3])
        print ens_id

        for l in open(f):
            fields = l.rstrip().split('\t')
            # The cluster contents
            cluster = eval(fields[3])
            print type(cluster[0]), cluster
            fields = fields[:3] + [ min(cluster), max(cluster) ] + fields[4:]
            print >>outfile, '\t'.join([ ens_id, dataset ] + [ str(field) for field in fields ])


if __name__ == "__main__":
    main()

