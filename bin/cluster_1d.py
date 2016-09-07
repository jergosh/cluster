import os
from os import path
import sys
import argparse
import pandas
import multiprocessing

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
    argparser.add_argument('--nthreads', metavar='nthreads', type=int, default=1)
    argparser.add_argument('--ret_thr', metavar='threshold', type=float, default=0.05)

    argparser.set_defaults(disrete=False)

    args = argparser.parse_args()

    
    slr = pandas.read_csv(open(args.infile), sep="\t", comment="\n")
    coords = [ [i] for i in range(0, slr.shape[0]) ]
    # find and remove NAs in slr?

    if args.discrete:
        marks = list((slr['Omega'] > 1) & (slr['Adj.Pval'] < args.thr))
        print marks
        if sum(marks) < 2:
            return
    else:
        marks = list(slr['Omega'])

    cluster_id = 1
    outfile = open(args.outfile, 'w')
    p = multiprocessing.Pool(args.nthreads)

    ids = range(len(marks))
    ret = cucala.signMWcont_multi(coords, marks, ids, cucala.order_dists(coords), args.niter, p)
    print >>outfile, '\t'.join([ str(it) for it in ret ] + [ str(cluster_id) ])

    while ret[4] < args.thr:
        cluster_id += 1

        # to_keep contains the indices of IDs which are not in the cluster
        # that has been found
        to_keep = [ i for i, item in enumerate(ids) if item not in ret[3] ]
        
        ids[:] = [ item for i, item in enumerate(ids) if i in to_keep ]
        coords[:] = [ item for i, item in enumerate(coords) if i in to_keep ]
        marks[:] = [ item for i, item in enumerate(marks) if i in to_keep ]

        ret = cucala.signMWcont_multi(coords, marks, ids, cucala.order_dists(coords), args.niter, pool)
        print >>outfile, '\t'.join([ str(it) for it in ret ] + [str(cluster_id)] )

    outfile.close()

if __name__ == "__main__":
    main()
