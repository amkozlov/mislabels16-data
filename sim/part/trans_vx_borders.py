#! /usr/bin/env python
import sys

# -------
# MAIN
# -------
if __name__ == "__main__":

    if len(sys.argv) < 4:
        print "Usage: ./trans_vx_borders.py AlignedEcoli.seq Ecoli.coord.txt AlignedEcoli.part.txt"
        sys.exit(0)

    offset = 0

    ecoli_align_fname = sys.argv[1]
    ecoli_part_fname = sys.argv[2]        
    output_fname = sys.argv[3]

    with open(ecoli_align_fname) as inf:
        align_seq = inf.readline().strip()
        
    ecoli_seq = align_seq.translate(None, "?.-")       
        
    print "Ecoli seq len unaligned / aligned: %d / %d\n" % (len (ecoli_seq), len (align_seq))
    
    ecoli_part = {}
    pos_map = {}
    with open(ecoli_part_fname) as inf:
        for line in inf:
            part_name, ranges = line.strip().split("=")
            part_ranges = []
            for range_str in ranges.split(","):
              range_min,range_max = [int(s)+offset for s in range_str.split("-")]
              part_ranges.append((range_min,range_max))
              pos_map[range_min] = None
              pos_map[range_max] = None
              
            ecoli_part[part_name] = part_ranges
            
    print "Original partitions:\n", ecoli_part

    ecoli_pos = 0
    align_pos = 0
    ecoli_start = None
    while ecoli_pos < len(ecoli_seq):
        while align_seq[align_pos] in ['?', '-', '.']:
           align_pos += 1
        if not ecoli_start:
            ecoli_start = align_pos
        if align_seq[align_pos] != ecoli_seq[ecoli_pos]:
            print "Error: base mismatch at alignment position %d: %s <> %s" % (align_pos+1, align_seq[align_pos], ecoli_seq[ecoli_pos])
            sys.exit()
        if ecoli_pos+1 in pos_map:
            pos_map[ecoli_pos+1] = align_pos+1
        ecoli_pos += 1
        align_pos += 1
        
    ecoli_end = align_pos+1
        
    print pos_map

    align_part = {}
    v_minmax = [ecoli_start, ecoli_end]
    for part_name, ecoli_ranges in ecoli_part.iteritems():
        align_ranges = []
        for r_min, r_max in ecoli_ranges:
            align_min = pos_map[r_min]
            align_max = pos_map[r_max]
            align_ranges.append((align_min, align_max))
            v_minmax += [align_min, align_max]
            
        align_part[part_name] = align_ranges
    
    cons_range = []
    lastmax_pos = None
    for pos in sorted(v_minmax):
        if lastmax_pos:
            cons_range += [(lastmax_pos+1, pos-1)]
            lastmax_pos = None
        else:
            lastmax_pos = pos

    align_part["cons"] = cons_range
    
    if ecoli_start > 1:
        align_part["flanks"] = [(1, ecoli_start), (ecoli_end, len(align_seq))]
    
    print "\nAligned partitions:\n", align_part

    
    with open(output_fname, "w") as outf:
        for part_name in sorted(align_part.keys()):
            outstr = "DNA, %s=%s\n" % (part_name, ",".join(["%d-%d" % (rmin, rmax) for rmin, rmax in align_part[part_name]]))
            outf.write(outstr);
