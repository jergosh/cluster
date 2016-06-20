import os
from os import path
import pandas
import sys
import argparse

import cucala
from slr import *


def main():
    argparser = argparse.ArgumentParser()

    argparser.add_argument('--infile', metavar='slr_file', type=str, required=True)
    argparser.add_argument('--niter', metavar='niter', type=int, default=999)

    args = argparser.parse_args()

    slr = pandas.read_csv(open(args.infile), sep="\t", comment="\n")
    coords = [ [i] for i in range(0, slr.shape[0]) ]
    # find and remove NAs in slr?

    marks = slr['Omega']

    print cucala.signMWcont(coords, marks, cucala.order_dists(coords), args.niter)

if __name__ == "__main__":
    main()
