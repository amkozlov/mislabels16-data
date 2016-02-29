#!/usr/bin/env python

import os
import sys
import shutil
import argparse
from operator import itemgetter
from epatax_common import EpataxConfig
from epatax_eval import EpataxEvaluator
from epatax_taxonomy import Taxonomy,GGTaxonomyFile


suffix_list = ["k__", "p__", "c__", "o__", "f__", "g__", "s__"]

def parseArgs():
    parser = argparse.ArgumentParser(description="Calculate accuracy statistics for mislabels output")
    parser.add_argument("-m", dest="mis_fname", required=True,
            help="""Mislabels file""")
    parser.add_argument("-t", dest="taxonomy_fname", default=None,
            help="""Full taxonomy with mislabeled taxa""")
    parser.add_argument("-g", dest="gtruth_fname", default=None,
            help="""Ground-truth taxonomy""")
    parser.add_argument("-o", dest="output_dir", default="",
            help="""Directory where output file should be created (default: directory of assignment file)""")
    parser.add_argument("-c", dest="config_fname", default="epatax.cfg",
            help="config file name (default: epatax.cfg)")
    parser.add_argument("-x", dest="min_conf", type=float, default=0.0,
            help="Confidence threshold (default: 0.0)")
    args = parser.parse_args()
    return args

def get_rank_level(rank_name):
        return { "Kingdom": 1,
                 "Phylum": 2,
                 "Class": 3,
                 "Order": 4,
                 "Family": 5,
                 "Genus": 6,
                 "Species": 7
                }[rank_name]

def get_rank_level_by_suffix(suffix):
        return { "k__": 1,
                 "p__": 2,
                 "c__": 3,
                 "o__": 4,
                 "f__": 5,
                 "g__": 6,
                 "s__": 7
                }[suffix]


def parse_mis_rec(line):
    toks = line.split("\t")
    if len(toks) < 2:
       toks = line.split()
    mis_rec = {}
    mis_rec["name"] = toks[0]
    mis_rec["mis_rank"] = toks[1]
    rank_suffix = toks[2][:3]
    if rank_suffix in suffix_list:
      mis_rec["mis_lvl"] = get_rank_level_by_suffix(rank_suffix)
    else:
      mis_rec["mis_lvl"] = get_rank_level(toks[1])
   
    mis_rec["conf"] = float(toks[4])
    mis_rec["orig_ranks"] = toks[5].split(";")
    mis_rec["ranks"] = toks[6].split(";")
    
    return mis_rec

def parse_mislabels(mis_fname, min_conf):
    mislabels = {}
    with open(mis_fname) as inf:
        last_rank = ""
        rank_lvl = 0
        for line in inf:
            if line.startswith(";"):
                continue
            mis_rec = parse_mis_rec(line)
            if mis_rec['mis_rank'] != last_rank:
                last_rank = mis_rec['mis_rank']
                rank_lvl += 1
            mis_rec['mis_rank'] = "%d_%s" % (rank_lvl, mis_rec['mis_rank'])
            # ignore species level
            if mis_rec['mis_lvl'] != 7 and mis_rec['conf'] >= min_conf:
                mislabels[mis_rec['name']] = mis_rec
    
    return mislabels

def read_assignment_file(assign_fname):
    tax_assign_map = {}
    tax_assign_conf = {}
    empty_ranks = GGTaxonomyFile.rank_placeholders + ["-","None", "No blast hit", "Unclassified"]
    with open(assign_fname) as fin:
        for line in fin:
            line = line.strip()
            toks = line.split("\t")
            
            # TODO fix unaligned sequences
            if len(toks) < 2:
                continue
                
            sid = toks[0]
            ranks_str = toks[1]
            ranks = ranks_str.split(";")
            for i in range(len(ranks)):
                rank_name = ranks[i].strip()
                if rank_name in empty_ranks:
                    rank_name = ""
                ranks[i] = rank_name
                
            if len(ranks) < 7:
                ranks += [""] * (7 - len(ranks))
            tax_assign_map[sid] = ranks     

            if len(toks) > 2 and ranks[0]:
                cf_toks = toks[2].split(";")
                cf_vals = []
                for cf_str in cf_toks:
                    cf_vals += [float(cf_str)]
                tax_assign_conf[sid] = cf_vals
            else:
                tax_assign_conf[sid] = [1.] * len(ranks)
                
    return tax_assign_map, tax_assign_conf


# -------
# MAIN
# -------
if __name__ == "__main__":
    args = parseArgs()
    config = EpataxConfig(args.config_fname, "")
#    if not os.path.isdir(config.results_dir):
#        os.mkdir(config.results_dir)

    if args.min_conf > 1.0:
      args.min_conf /= 100.

    print "Mislabels with confidence >= %f:" % args.min_conf 

#    (tax_assign_map, tax_assign_conf) = read_assignment_file(args.assign_fname)
    mislabel_map = parse_mislabels(args.mis_fname, args.min_conf)

    rank_mis_map = {}
    for mis_rec in mislabel_map.itervalues():
        rank = mis_rec["mis_rank"]
        rank_mis_map[rank] = rank_mis_map.get(rank, 0) + 1

    rank_mis_count = sorted(rank_mis_map.items(), key=itemgetter(0))

    total = 0
    for rank, cnt in rank_mis_count:
        total += cnt
#        print "%s:\t%d\t%d" % (rank,cnt,total)
    print ""

    if not args.gtruth_fname:
        sys.exit(0)
    
    true_taxonomy = GGTaxonomyFile(args.gtruth_fname)
    
    out_path, out_stem = os.path.split(args.mis_fname)
    if args.output_dir:
        out_path = args.output_dir
    if out_stem.endswith(".mis"):
        out_stem = out_stem[:-4]
#    mis_fname = os.path.join(out_path, out_stem + ".mis")
    
    e = EpataxEvaluator(config, args.taxonomy_fname)
    e.calc_mislabel_stats(mislabel_map, true_taxonomy.seq_ranks_map, False)
    
    fp_fname = os.path.join(out_path, out_stem + ".eval%d.fp" % int(args.min_conf * 100))
    with open(fp_fname, "w") as fout:
        for sid in e.fp:
            ranks = mislabel_map[sid]["ranks"]
            lvl = mislabel_map[sid]["mis_rank"]
            true_ranks = e.taxonomy.get_seq_ranks(sid)
            fout.write("%s\t%s\t%s\t%s\n" % (sid, lvl, Taxonomy.lineage_str(true_ranks), Taxonomy.lineage_str(ranks)))

    tp_fname = os.path.join(out_path, out_stem + ".eval%d.tp" % int(args.min_conf * 100))
    with open(tp_fname, "w") as fout:
        for sid in e.tp:
            ranks = mislabel_map[sid]["ranks"]
            fout.write("%s\t%s\n" % (sid, Taxonomy.lineage_str(ranks)))

    fn_fname = os.path.join(out_path, out_stem + ".eval%d.fn" % int(args.min_conf * 100))
    with open(fn_fname, "w") as fout:
        for sid in e.fn:
            true_ranks = true_taxonomy.get_seq_ranks(sid)
            ranks = e.taxonomy.get_seq_ranks(sid)
            fout.write("%s\t%s\t%s\n" % (sid, Taxonomy.lineage_str(true_ranks), Taxonomy.lineage_str(ranks)))

    full_fname = os.path.join(out_path, out_stem + ".eval%d.full" % int(args.min_conf * 100))
    with open(full_fname, "w") as fout:
#        for sid in mislabel_map:
#            ranks = mislabel_map[sid]["ranks"]
#            fout.write("%s\t%s\n" % (sid, Taxonomy.lineage_str(ranks)))
        for sid in e.tp:
            ranks = mislabel_map[sid]["ranks"]
            fout.write("%s\t%s\n" % (sid, Taxonomy.lineage_str(ranks)))
        for sid in e.fn:
            ranks = e.taxonomy.get_seq_ranks(sid)
            fout.write("%s\t%s\n" % (sid, Taxonomy.lineage_str(ranks)))

    stats_fname = os.path.join(out_path, out_stem + ".eval%d.ident" % int(args.min_conf * 100))
    e.write_stats(stats_fname)
#    e.write_misclassified(tax_assign_map, tax_assign_conf, mis_fname)
    print "Statistics has been written to: %s" % stats_fname
