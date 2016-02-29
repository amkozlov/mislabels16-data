#!/usr/bin/env python
import sys

if __name__ == "__main__":
    if len(sys.argv) < 3: 
        print "Usage: ./merge.py ali1.phy ali2.phy ..."
        sys.exit()

    total_sites = 0
    out_phy = {}

    for phy_fname in sys.argv[1:]:
      with open(phy_fname) as inf:
        header = inf.readline()
        (ntaxa, nsites) = header.split()
        nsites = int(nsites)

        total_sites += nsites

        gaps = [0] * nsites

        while True:
            line = inf.readline()
            if not line or not line.strip():
              break
            (taxon, seq) = line.strip().split()
            if taxon in out_phy:
              out_phy[taxon] += seq.strip()
            else:
              out_phy[taxon] = seq.strip()


    print "%d %d" % (len(out_phy), total_sites)
    for taxon, seq in out_phy.iteritems():
      print "%s %s" % (taxon, seq)