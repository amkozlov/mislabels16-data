#! /usr/bin/env python

import sys

if __name__ == "__main__":
    if len(sys.argv) == 1: 
        print "Usage: filter_tax.py name2acc.txt id2acc.txt fulltax.txt"
        sys.exit()

    map_fname = sys.argv[1]
    gg_map_fname = sys.argv[2]
    tax_fname = sys.argv[3]
    name2acc = {}
    acc2name = {}
    with open(map_fname) as mapf:
	for line in mapf:
            name, acc = line.split()
	    name2acc[name] = acc
	    acc2name[acc] = name

    ggid2acc = {}
    acc2ggid = {}
    with open(gg_map_fname) as mapf:
        for line in mapf:
            ggid, acc = line.split()
            acc = acc.split(".", 1)[0]
            ggid2acc[ggid] = acc
            acc2ggid[acc] = ggid

    acc_list = set(acc2name.keys())
    with open(tax_fname) as inf:
	for line in inf:
	    ggid, rest = line.strip().split(None, 1)
            acc = ggid2acc.get(ggid, None)
	    if acc in acc_list:
		name = acc2name.get(acc, acc)
		toks = rest.strip().split(None,1)
		print "%s\t%s" % (name, rest)
		if acc in acc_list:
		    acc_list.remove(acc)

#    if len(acc_list) > 0:
#      print "\nNot found: %d\n" % len(acc_list)
#      for acc in acc_list:
#	print "%s\t%s" % (acc2name.get(acc, acc), "Unassigned")
