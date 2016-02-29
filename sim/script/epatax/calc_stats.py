#!/usr/bin/env python

import os
import sys
import shutil
import argparse
from epatax_common import EpataxConfig
from epatax_eval import EpataxEvaluator
from epatax_taxonomy import GGTaxonomyFile

def parseArgs():
    parser = argparse.ArgumentParser(description="Calculate assignment accuracy statistics")
    parser.add_argument("-a", dest="assign_fname", required=True,
            help="""Taxonomic assignment file, i.e. output of EPA classifier or other taxonomic placement program""")
    parser.add_argument("-t", dest="taxonomy_fname", required=True,
            help="""Ground-truth taxonomy""")
    parser.add_argument("-o", dest="output_dir", default="",
            help="""Directory where output file should be created (default: directory of assignment file)""")
    parser.add_argument("-c", dest="config_fname", default="epatax.cfg",
            help="config file name (default: epatax.cfg)")
    args = parser.parse_args()
    return args


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

    (tax_assign_map, tax_assign_conf) = read_assignment_file(args.assign_fname)
    
    out_path, out_stem = os.path.split(args.assign_fname)
    if args.output_dir:
        out_path = args.output_dir
    if out_stem.endswith(".full"):
        out_stem = out_stem[:-5]
    if out_stem.endswith(".tp") or out_stem.endswith(".fp"):
        out_stem = out_stem[:-3]
    stats_fname = os.path.join(out_path, out_stem + ".corr")
    mis_fname = os.path.join(out_path, out_stem + ".mis")
    
    e = EpataxEvaluator(config, args.taxonomy_fname)
    e.calc_stats(tax_assign_map, False)
    e.write_stats(stats_fname)
    e.write_misclassified(tax_assign_map, tax_assign_conf, mis_fname)
    print "Statistics has been written to: %s" % stats_fname
