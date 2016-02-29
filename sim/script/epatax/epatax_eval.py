#!/usr/bin/env python

from epatax_taxonomy import GGTaxonomyFile, Taxonomy

class StatEntry:
    
    def __init__(self): 
        self.tp = 0
        self.fp = 0
        self.fp2 = 0
        self.tn = 0
        self.fn = 0
        self.precision = 0
        self.precision2 = 0
        self.recall = 0
        self.f1 = 0

    def recalc(self):
        if self.tp + self.fn == 0 or self.tp + self.fp == 0:
            return
        self.precision = float(self.tp) / (self.tp + self.fp)
        self.precision2 = float(self.tp) / (self.tp + self.fp + self.fp2)
        self.recall = float(self.tp) / (self.tp + self.fn)
        if self.precision + self.recall > 0:
            self.f1 = 2 * self.precision * self.recall / (self.precision + self.recall)

    @staticmethod    
    def desc_string():
#        return "\t".join([";Level", "TP", "FP", "FP2", "TN", "FN", "PRE", "PRE2", "REC", "F1"])
         return "\t".join([";Level", "TP", "FP", "TN", "FN", "PRE", "REC", "F1"])
    
    def to_string(self, delim="\t"):
#        cnt_vals = [self.tp, self.fp, self.fp2, self.tn, self.fn]
#        acc_vals = [self.precision, self.precision2, self.recall, self.f1]
        cnt_vals = [self.tp, self.fp, self.tn, self.fn]
        acc_vals = [self.precision, self.recall, self.f1]
        return delim.join(map(str, cnt_vals) + map("{:5.3f}".format, acc_vals))

    def from_string(self, s, delim="\t"):
        vals = s.split(delim)
        self.tp = int(vals[0])
        self.fp = int(vals[1])
        self.tn = int(vals[2])
        self.fn = int(vals[3])

    def add(self, other):
        self.tp += other.tp
        self.fp += other.fp
        self.fp2 += other.fp2
        self.fn += other.fn
        self.tn += other.tn
    
class EpataxEvaluator:
    
    def __init__(self, config, taxonomy_fname): 
        self.config = config
        self.taxonomy = GGTaxonomyFile(taxonomy_fname)
        self.max_level = self.taxonomy.max_rank_level()
        self.stats = [ StatEntry() for i in range(self.max_level+1) ]
        self.totstat_lvl = None
        self.miss_list = []
        self.tp = []
        self.tn = []
        self.fp = []
        self.fn = []
    
    def set_taxonomy(self, taxonomy_fname):
        self.taxonomy = GGTaxonomyFile(taxonomy_fname)

    def calc_stats(self, tax_assign_map, cumulative=True):
        if not cumulative:
          self.stats += [ StatEntry() ]
          self.totstat_lvl = len(self.stats)-1

        t = self.taxonomy        
        # evaluate the classification results against the ground truth
        for sid, epa_ranks in tax_assign_map.items():
            true_ranks = t.get_seq_ranks(sid)
            r = self.max_level
            miss = False     

            # for clades which were ingored during the reference tree construction,
            # true assignment will be the empty one (i.e., if we excluded a genus G1
            # from the tree, the "right" assignment will be the family it belongs to)
            ignore_level = len(true_ranks)
            for (rank_level, rank_name) in self.config.eval_ignored_clades:
                if rank_level < ignore_level and true_ranks[rank_level] == rank_name:
                    ignore_level = rank_level

            while r >= 0:

                # as described above, we treat ranks from ignored clades as 
                # "missing data" in the ground truth
                if r >= ignore_level:
                    # NOTE: even though no assignment is a technically TN, we consider 
                    # it being a TP in order to calculate precision/recall as usual
                    if epa_ranks[r] == "":
                        self.stats[r].tp += 1
                    else:
                        self.stats[r].fp += 1
                # "normal" clades
                elif epa_ranks[r] == "" and true_ranks[r] == "":
                    self.stats[r].tn += 1
                elif epa_ranks[r] == "" and true_ranks[r] != "":
                    self.stats[r].fn += 1
                    miss = True
                elif epa_ranks[r] != "" and true_ranks[r] == "":
                    # TODO think about this case
                    # it's considered FP for now, but can be actually correct classification 
                    self.stats[r].fp2 += 1
                    miss = True
                elif epa_ranks[r] == true_ranks[r]:
                    self.stats[r].tp += 1
                elif epa_ranks[r] != true_ranks[r]:
                    self.stats[r].fp += 1
                    miss = True
                else:
                    print "FATAL ERROR:  Oops, this shouldn't have happened..."
                    sys.exit()
                r -= 1
            if miss:
                self.miss_list += [sid]
        
        # calculate accuracy measures
        for i in range(len(self.stats)):
          self.stats[i].recalc()          

    def calc_ranktest_stats(self, tax_assign_map):
        t = self.taxonomy        
        # evaluate the classification results against the ground truth
        for sid, epa_ranks in tax_assign_map.items():
            true_ranks = t.get_seq_ranks(sid)
            r = self.max_level
            miss = False     
            fp = False

            # for clades which were ingored during the reference tree construction,
            # true assignment will be the empty one (i.e., if we excluded a genus G1
            # from the tree, the "right" assignment will be the family it belongs to)
            ignore_level = len(true_ranks)            
            for (rank_level, rank_name) in self.config.eval_ignored_clades:
                if rank_level < ignore_level and true_ranks[rank_level] == rank_name:
                    ignore_level = rank_level

            while r >= 0:

                # as described above, we treat ranks from ignored clades as 
                # "missing data" in the ground truth
                if r >= ignore_level:
                    # NOTE: even though no assignment is a technically TN, we consider 
                    # it being a TP in order to calculate precision/recall as usual
                    if epa_ranks[r] == "":
                        self.stats[r].tp += 1
                    else:
                        self.stats[r].fp += 1
                        fp = True
                # FP in lower rank -> FP in all upper ranks
                elif fp:
                    self.stats[r].fp += 1
                # "normal" clades
                elif epa_ranks[r] == "" and true_ranks[r] == "":                
                    self.stats[r].tn += 1
                elif epa_ranks[r] == "" and true_ranks[r] != "": 
                    self.stats[r].fn += 1
                    miss = True
                elif epa_ranks[r] != true_ranks[r]:
                    self.stats[r].fp += 1
                    miss = True
                    fp = True
                elif epa_ranks[r] != "" and true_ranks[r] == "": 
                    # TODO think about this case
                    # it's considered FP for now, but can be actually correct classification 
                    self.stats[r].fp2 += 1
                    miss = True
                    fp = True
                elif epa_ranks[r] == true_ranks[r]:
                    self.stats[r].tp += 1
                else:
                    print "FATAL ERROR:  Oops, this shouldn't have happened..."
                    sys.exit()
                r -= 1
            if miss:
                self.miss_list += [sid]
        
        # calculate accuracy measures
        for i in range(len(self.stats)):
          self.stats[i].recalc()   
          
    def get_mis_level(self, ranks, true_ranks):
        mislabel_lvl = -1
        min_len = min(len(ranks),len(true_ranks))
        for rank_lvl in range(min_len):
            if ranks[rank_lvl] != "-" and ranks[rank_lvl] != true_ranks[rank_lvl]:
                mislabel_lvl = rank_lvl
                break
        
        return mislabel_lvl
    
    def calc_mislabel_stats(self, mislabel_map, true_map, cumulative=True):
        t = self.taxonomy        

        if not cumulative:
          self.stats += [ StatEntry() ]
          self.totstat_lvl = len(self.stats)-1
        
        for sid in mislabel_map.iterkeys():
            if sid not in true_map:
                self.fp += [sid]
                mis_rec = mislabel_map[sid]
                lvl = mis_rec["mis_lvl"]
                rlist = range(lvl-1, 7) if cumulative else [lvl-1, self.totstat_lvl]
                for r in rlist:
                    self.stats[r].fp += 1
            
        # evaluate the classification results against the ground truth
        for sid in true_map.iterkeys():
            true_ranks = true_map[sid]
            lvl = self.get_mis_level(t.get_seq_ranks(sid), true_ranks)
            rlist = range(lvl, 7) if cumulative else [lvl, self.totstat_lvl]
            if sid not in mislabel_map:
                self.fn += [sid]
                for r in rlist:
                    self.stats[r].fn += 1
            else:
                self.tp += [sid]
                mis_rec = mislabel_map[sid]
#                lvl = mis_rec["mis_lvl"]
                for r in rlist:
                    self.stats[r].tp += 1

        # calculate accuracy measures
        for i in range(len(self.stats)):
          self.stats[i].recalc()

    def write_misclassified(self, tax_assign_map, tax_assign_conf, out_fname):
        t = self.taxonomy        
        with open(out_fname, 'w') as fout:
            fout.write("; List of misclassified sequences.\n")
            fout.write("; Format: 4 lines per sequence\n")
            fout.write("; 1 sequence ID\n")
            fout.write("; 2 correct assignment\n")
            fout.write("; 3 epatax assignment\n")
            fout.write("; 4 rank confidence levels\n\n")

            for sid in self.miss_list:
                fout.write(sid + "\n")
                true_ranks = t.get_seq_ranks(sid)
                fout.write(Taxonomy.lineage_str(true_ranks) + "\n")
                epa_ranks = tax_assign_map[sid]
                fout.write(Taxonomy.lineage_str(epa_ranks) + "\n")
                rank_conf = tax_assign_conf[sid]
                fout.write("\t".join(["%.3f" % conf for conf in rank_conf]) + "\n")
                fout.write("\n")

    def write_stats(self, out_fname):
        t = self.taxonomy        
        with open(out_fname, 'w') as fout:
            fout.write(StatEntry.desc_string() + "\n")
            for i in range(len(self.stats)):
               rname = t.rank_level_name(i) if i != self.totstat_lvl else "Total"
               fout.write(rname + "\t" + self.stats[i].to_string("\t") + "\n")
