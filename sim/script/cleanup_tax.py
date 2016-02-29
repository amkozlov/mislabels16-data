#!/usr/bin/env python
import sys

if __name__ == "__main__":
    if len(sys.argv) == 1: 
        print "Usage: ./cleanup_tax.py input.tax"
        sys.exit()

    intax_fname = sys.argv[1]

    prefix = ["k__", "p__", "c__", "o__", "f__", "g__", "s__"]

    with open(intax_fname, 'r') as inf:
      for line in inf:
        (seqid, taxpath) = line.split("\t")
        taxranks = taxpath.strip().split(";")
        newranks = [""] * 7
        rankno = 0
        for rank in taxranks[:-2]:
          if len(taxranks) > 7 and rank.endswith(("ineae", "idae")):
            continue
          newranks[rankno] = prefix[rankno] + rank.replace(" ", "_")
          rankno += 1
          if rankno == 5:
            break

        incomplete = (rankno < 5)

        while rankno < 5:
          newranks[rankno] = prefix[rankno] + taxranks[-2].replace(" ", "_")
          rankno += 1

        newranks[-2] = prefix[-2] + taxranks[-2].replace(" ", "_")
        newranks[-1] = prefix[-1] + taxranks[-1].replace(" ", "_")
        rankno += 2

        if incomplete:
          print "%s\t%s\t%s" % (seqid, taxpath.strip(), ";".join(newranks))

        if rankno != 7:
            print "Problem!"
            print taxranks, newranks
            sys.exit()

#        print "%s\t%s" % (seqid, ";".join(newranks))
