#! /usr/bin/env python

import sys

if __name__ == "__main__":
    if len(sys.argv) == 1: 
        print "Usage: acc2name.py src_file.txt name2acc_map.txt [-r]"
        sys.exit()

    in_fname = sys.argv[1]
    map_fname = sys.argv[2]
    if len(sys.argv) == 4 and sys.argv[3] == "-r":
	reverse = True
    else:
	reverse = False

    name2acc = {}
    acc2name = {}
    acc2name_s = {}
    with open(map_fname) as mapf:
	for line in mapf:
	    name, acc = [i.strip() for i in line.split()]
#            print "|%s|%s|" % (name, acc)
	    name2acc[name] = acc
	    acc2name[acc] = name
            if acc[-2] == '.':
              acc2name_s[acc[:-2]] = name

    with open(in_fname) as inf:
	for line in inf:
	    acc, rest = line.split(None, 1)
 #           print "|%s|" % acc
	    if reverse:
		name = name2acc.get(acc, None)
	    else:
		name = acc2name.get(acc, None)
                if not name:
                  name = acc2name_s.get(acc, None)
            if not name:
              name = acc
	    print "%s\t%s" % (name, rest.strip())
