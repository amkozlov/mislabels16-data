# randomly introduce mislabels into taxonomy
# by J.Zhang, modified by A.Kozlov, HITS gGmbH

#!/usr/bin/env python
import os
import sys
import random
import math
import copy

class rd_generator:
    def __init__(self, seed):
        self.rand_nr = random.Random()
        self.rand_nr.seed(seed)

def map2string(m):
    outs = ""
    for key in m:
        v = m[key]
        #print(key)
        #print(v)
        vv = ""
        for vi in v:
            vv = vv + vi + ";"
        vv = vv[0:-1]
        outs = outs + key + "	" + vv + "\n"
    return outs

def rank2string(l):
    s = ""
    for e in l:
        s = s + e + ";"
    s = s[0:-1]
    return s

def is_same_ranks(rk1, rk2, k):
    flag = True
    for i in range(k):
        if rk1[i] != rk2[i]:
            flag = False
    return flag


class rd_mislable:
    def __init__(self, tax_input, seed, p_mis = 0.01):
        self.tax = {}
        self.names = []
        prefix = ["k__", "p__", "c__", "o__", "f__", "g__", "s__"]
        with open(tax_input) as fin:
            for line in fin:
                ll = line.strip().split("\t")
                taxonomy = ll[1].split(";")
                if len(taxonomy) < 7:
                  taxnew = [""] * 7
                  taxnew[6] = taxonomy[-1]
                  taxnew[5] = taxonomy[-2]
                  for i in range(5):
                    if i < len(taxonomy)-2:
                      taxnew[i] = taxonomy[i]
                    else:
                      taxnew[i] = prefix[i] + taxonomy[-2]
                  print ll[0], "\t", taxnew
                  self.tax[ll[0]] = taxnew
                else:
                  self.tax[ll[0]] = taxonomy
                self.names.append(ll[0])
        self.num_taxa = len(self.names)
        self.num_mis = int(self.num_taxa * p_mis)
       #self.prob = [0.05, 0.05, 0.1, 0.3, 0.5] #p, c, o, f ,g.
        self.prob = [0.05, 0.10, 0.15, 0.35, 0.35] #p, c, o, f ,g.
#        self.prob = [0.00, 0.1, 0.2, 0.2, 0.5] #p, c, o, f ,g.
        self.num = [int(x*self.num_mis) for x in self.prob]
        print "Total taxa: ",self.num_taxa, ", mislabels by rank: ",  self.num
        self.rd = rd_generator(seed = seed)
        self.truetax = {}
        self.mistax = {}
        self.remaintax = {}
        self.testingtax = {}
        
    def find_all_taxs(self, rank_idx, rank_name_exclude, ranks2keep):
        rks = []
        for key in self.tax:
            rank = self.tax[key]
            if rank[rank_idx] != self.tax[rank_name_exclude][rank_idx]:
                if is_same_ranks(ranks2keep, rank, rank_idx):
                    #print(rank[rank_idx])
                    #print("vs")
                    #print(self.tax[rank_name_exclude][rank_idx])
                    rks.append(rank)
        return rks
    
    def mis_lable(self, names, rank_idx):
        for name in names:
            remain_ranks = self.find_all_taxs(rank_idx, name, self.tax[name])
            if len(remain_ranks) > 0:
                idx = self.rd.rand_nr.randint(0, len(remain_ranks)-1)
                tk1 = self.tax[name]
                tk2 = copy.copy(tk1)
                tk2.append(str(rank_idx))
                self.truetax[name] = tk2
                self.mistax[name] = remain_ranks[idx]
            else:
                self.remaintax[name] = self.tax[name]
                print "find single rank: ", self.tax[name][rank_idx]
    
    def simulate(self, fout):
        idx = range(len(self.names))
        self.rd.rand_nr.shuffle(idx)
        names1 = []
        names2 = []
        names3 = []
        names4 = []
        names5 = []
        names_unchanged = []
        names_rest = []
        cnt = 0 
        for i in range(self.num[0]):
            names1.append(self.names[idx[cnt]])
            cnt = cnt + 1
        
        for i in range(self.num[1]):
            names2.append(self.names[idx[cnt]])
            cnt = cnt + 1        

        for i in range(self.num[2]):
            names3.append(self.names[idx[cnt]])
            cnt = cnt + 1
        
        for i in range(self.num[3]):
            names4.append(self.names[idx[cnt]])
            cnt = cnt + 1

        for i in range(self.num[4]):
            names5.append(self.names[idx[cnt]])
            cnt = cnt + 1
            
        for i in range(self.num_mis):
            names_unchanged.append(self.names[idx[cnt]])
            cnt = cnt + 1
        
        while cnt < len(self.names):
            names_rest.append(self.names[idx[cnt]])
            cnt = cnt + 1

        self.mis_lable(names1, rank_idx = 1)
        self.mis_lable(names2, rank_idx = 2)
        self.mis_lable(names3, rank_idx = 3)
        self.mis_lable(names4, rank_idx = 4)
        self.mis_lable(names5, rank_idx = 5)
        
        for name in names_unchanged:
            self.testingtax[name] = self.tax[name]
        
        for name in names_rest:
            self.remaintax[name] = self.tax[name]
            
        s_mis = map2string(self.mistax)
        s_mis_true = map2string(self.truetax)
        s_testing = map2string(self.testingtax)
        s_remain = map2string(self.remaintax)
        
        with open(fout+".tax", "w") as fo:
            fo.write(s_mis + s_testing + s_remain)
        
        #with open(fout+".testing.tax", "w") as fo:
            #fo.write(s_mis)
            
        with open(fout+".true.tax", "w") as fo:
            fo.write(s_mis_true)
        

if __name__ == "__main__":
    
    if len(sys.argv) < 2:
      print "usage: python tax_simulator.py <simdir> <pct> <repnum>"
      sys.exit()

    sim_dir = sys.argv[1] + "/"
    pct = int(sys.argv[2])
    repnum = int(sys.argv[3])

    input_tax = sim_dir + "p0/true.tax"
    out_dir = sim_dir + "p%d" % pct
    out_file = out_dir + "/simfull-p%d-t%d" % (pct, repnum)
     
    rm = rd_mislable(input_tax, seed = repnum, p_mis = 0.01 * pct)
    rm.simulate(fout = out_file)
   
    #rm.simulate(fout = "/home/zhangje/GIT/tax_benchmark/simulator/mLTP10")
    
