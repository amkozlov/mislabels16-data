#! /usr/bin/env python
import sys
import os

from epatax_eval import StatEntry

if __name__ == "__main__":
    if len(sys.argv) == 1: 
        print "Usage: summarize.py PCT sativa|rdp|uclust ident|corr CUTOFF"
        sys.exit()
    
    pct = int(sys.argv[1])
    method = sys.argv[2]
    metric = sys.argv[3]
    cutoff = int(sys.argv[4])
    
    ranks = ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Total"]
    total_stats = []
    for r in ranks:
        total_stats.append(StatEntry())
    
    rep = 1
    while os.path.isdir("t%d" % rep):
        stats_fname = "t%d/%s-p%d-t%d.eval%d.%s" % (rep, method, pct, rep, cutoff, metric)
        rank_num = 0
        with open(stats_fname) as inf:
            for line in inf:
                if not line.startswith(";"):
                    rank, stats = line.split("\t", 1)
                    entry = StatEntry()
                    entry.from_string(stats)
                    total_stats[rank_num].add(entry)
                    rank_num += 1
        rep += 1

    out_fname = "%s-p%d-total.eval%d.%s" % (method, pct, cutoff, metric)
    with open(out_fname, "w") as outf:
        outf.write(StatEntry.desc_string() + "\n")
        for i in range(len(ranks)):
            total_stats[i].recalc()
            outs = total_stats[i].to_string()
            outf.write(ranks[i] + "\t" + outs + "\n")

    print "Processed %d files, output is saved to: %s\n" % (rep-1, out_fname)
