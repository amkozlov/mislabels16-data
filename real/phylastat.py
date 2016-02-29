#! /usr/bin/env python

import sys
import operator

if __name__ == "__main__":
    if len(sys.argv) == 1: 
        print "Usage: phylastat.py full.tax input.mis [cutoff]"
        sys.exit()

    tax_fname = sys.argv[1]
    mis_fname = sys.argv[2]

    if len(sys.argv) > 3:
	cutoff = float(sys.argv[3])
    else:
	cutoff = 0.

    main_phyla = ["Proteobacteria", "Firmicutes", "Bacteroidetes", "Actinobacteria", "p__Proteobacteria", "p__Firmicutes", "p__Bacteroidetes", "p__Actinobacteria"]

    OTHERS="__Other"

    t_count = {}
    with open(tax_fname) as taxf:
	for line in taxf:
	    if line[0] == ";":
		continue
	    toks = line.split("\t")
	    orig_ranks = toks[1].split(";")
	    if len(orig_ranks) > 1:
		phylum = orig_ranks[1].strip().strip('"')
		t_count[phylum] = t_count.get(phylum, 0) + 1
		if not phylum in main_phyla:
		    t_count[OTHERS] = t_count.get(OTHERS, 0) + 1
	    else:
		print "No phylum in ranks: %s" % ";".join(orig_ranks)

    t_sorted = sorted(t_count.items(), key=operator.itemgetter(0))

    for name, cnt in t_sorted:
	print "%s\t%d" % (name, cnt)
    print "\n"

    m_count = {}
    with open(mis_fname) as mapf:
	for line in mapf:
	    if line[0] == ";":
		continue
	    toks = line.split("\t")
	    conf = float(toks[4])
	    if conf < cutoff:
		continue
	    orig_ranks = toks[5].split(";")
	    if len(orig_ranks) > 1:
		phylum = orig_ranks[1].strip().strip('"')
		m_count[phylum] = m_count.get(phylum, 0) + 1
		if not phylum in main_phyla:
		    m_count[OTHERS] = m_count.get(OTHERS, 0) + 1
	    else:
		print "No phylum in ranks: %s" % ";".join(orig_ranks)

    m_sorted = sorted(m_count.items(), key=operator.itemgetter(1), reverse=True)

    for name, cnt in m_sorted:
	if name in t_count:
	    pct = str(float(cnt) / t_count[name])
	else:
	    pct = "NA"
	print "%s\t%d\t%s" % (name, cnt, pct)