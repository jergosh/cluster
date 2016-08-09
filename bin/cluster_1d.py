import os
from os import path
import sys
import argparse
import pandas

import cucala
# Remove ?
from slr import *


def main():
    argparser = argparse.ArgumentParser()

    argparser.add_argument('--infile', metavar='slr_file', type=str, required=True)
    argparser.add_argument('--outfile', metavar='out_file', type=str, required=True)

    argparser.add_argument('--discrete', dest='discrete', action='store_true')
    argparser.add_argument('--continuous', dest='discrete', action='store_false')
    argparser.add_argument('--thr', metavar='threshold', type=float, default=0.05)
    argparser.add_argument('--niter', metavar='niter', type=int, default=999)
    argparser.add_argument('--ret_thr', metavar='threshold', type=float, default=0.05)
    argparser.set_defaults(disrete=False)

    args = argparser.parse_args()

    slr = pandas.read_csv(open(args.infile), sep="\t", comment="\n")
    coords = [ [i] for i in range(0, slr.shape[0]) ]
    # find and remove NAs in slr?

    if args.discrete:
        marks = list((slr['Omega'] > 1) & (slr['Adj.Pval'] < args.thr))
    else:
        marks = list(slr['Omega'])

    cluster_id = 1
    outfile = open(args.outfile, 'w')
    ret = cucala.signMWcont(coords, marks, cucala.order_dists(coords), args.niter)
    print >>outfile, '\t'.join([ str(it) for it in ret ] + [ str(cluster_id) ])

    while ret[4] < args.thr:
        cluster_id += 1
        coords[:] = [ item for i, item in enumerate(coords) if i not in ret[1] ]
        marks[:] = [ item for i, item in enumerate(marks) if i not in ret[1] ]

        ret = cucala.signMWcont(coords, marks, cucala.order_dists(coords), args.niter)
        print >>outfile, '\t'.join([ str(it) for it in ret ] + [str(cluster_id)] )

    outfile.close()

if __name__ == "__main__":
    main()
