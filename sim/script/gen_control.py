#!/usr/bin/env python

import sys

info_fname = sys.argv[1]

res_fname = info_fname.replace("info", "result", 1) + ".PARTITION"

i = 0

if len(sys.argv) > 2:
  p = sys.argv[2]
else:
  p = None

# 882	????????	30	105	64	106	57	57	56	51	12

# 45.52828	92.152951	75.820158	74.710209	47.626688	54.275103	51.135072	53.198551	42.909645	33.271214	75.165148

# 37.065931	87.653251	86.253521	34.680973	69.771679	42.990239	57.826729	65.244343	37.763092	39.211198	96.222209


ins_a = [ 2.0, 2.0 ] + [2.0] * 9
ins_M = [ 120, 100] + [15, 80, 25, 65, 35, 35, 28, 24, 7]
ins_rate = [ 0.0019, 0.0005 ] + [0.00175, 0.002, 0.0008, 0.0011, 0.001, 0.002, 0.0022, 0.0017, 0.0015]

root_len = [ 890, 45 ] + [ 28, 105, 55, 106, 57, 57, 56, 51, 15]

with open(info_fname, "r") as inf:
  while True:
    line = inf.readline()
    if not line:
      break;

    if not line.startswith("Model Parameters of Partition"):
      continue;

    pname = line.split(",")[1].split(": ")[1]
#    print pname

    if p and pname != p:
      i += 1
      continue

    alpha = float(inf.readline().split(": ")[1])
#    print alpha

    treelne = inf.readline().split(": ")[1]

    rAC = float(inf.readline().split(": ")[1])
    rAG = float(inf.readline().split(": ")[1])
    rAT = float(inf.readline().split(": ")[1])
    rCG = float(inf.readline().split(": ")[1])
    rCT = float(inf.readline().split(": ")[1])
    rGT = float(inf.readline().split(": ")[1])

    # re-normalize by rAG = f = 1
    rAC, rAG, rAT, rCG, rCT, rGT = ( r / rAG for r in (rAC, rAG, rAT, rCG, rCT, rGT) )

#    print rAC, rAG, rAT, rCG, rCT, rGT

    inf.readline()

    piA = float(inf.readline().split(": ")[1])
    piC = float(inf.readline().split(": ")[1])
    piG = float(inf.readline().split(": ")[1])
    piT = float(inf.readline().split(": ")[1])

    tree_fname = "%s.%d" % (res_fname, i)
    with open (tree_fname, "r") as treef:
      treestr = treef.read().strip()

    if treestr.endswith(":0.0;"):
      treestr = treestr[:-5] + ";"

#    print piA, piC, piG, piT

    print "[MODEL]    m_%s" % pname
    print "\t[submodel]  GTR %f %f %f %f %f" % (rCT, rAT, rGT, rAC, rCG)
    print "\t[statefreq] %f %f %f %f" % (piT, piC, piA, piG)
    print "\t[rates] 0 %f 0" % alpha


    print "\t[indelmodel]   LAV  %f %d" % (ins_a[i], ins_M[i])
    print "\t[indelrate]    %f" % ins_rate[i]

#    print "\t[insertmodel]   LAV  %f %d" % (ins_a[i], ins_M[i])
#    print "\t[insertrate]    %f" % ins_rate[i]
#    print "\t[deletemodel]   LAV  %f %f" % (ins_a[i], ins_M[i])
#    print "\t[deleterate]    0.0005"

    print ""
    print "[TREE] t_%s  %s" % (pname, treestr)
    print ""

    print "[PARTITIONS] p_%s   [t_%s m_%s %d]" % (pname, pname, pname, root_len[i])
    print ""

    print "[EVOLVE] p_%s 1 sim_%s" % (pname, pname)

    i += 1
