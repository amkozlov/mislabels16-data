#!/usr/bin/env python
import sys
import random
import numpy

if __name__ == "__main__":
    if len(sys.argv) == 1: 
        print "Usage: ./gaps.py ali.phy "
        sys.exit()

    phy_fname = sys.argv[1]

    nsites = None
    
    with open(phy_fname) as inf:
        header = inf.readline()
        (ntaxa, nsites) = header.split()
        nsites = int(nsites)

        gaps = [0] * nsites

        while True:
            line = inf.readline()
            if not line or not line.strip():
              break
            (taxon, seq) = line.strip().split()
            for i in range(len(seq)):
               if seq[i] in ["-", "?", "N"]:
                  gaps[i] += 1

    pgaps = [0] * nsites
    for i in range(nsites):
      pgaps[i] = float(gaps[i]) / float(ntaxa)

    hist, bin_edges = numpy.histogram(pgaps, range=(0.0, 1.0))
#    print hist, bin_edges
             
    total_base = float(ntaxa) * float(nsites)
    print "sites: %d, gappiness: %f %%" % (nsites, sum(gaps) / total_base * 100)
