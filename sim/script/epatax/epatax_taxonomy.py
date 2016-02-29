#!/usr/bin/env python

from ete2 import Tree
#from epatax_align import AlignmentUtils

class TaxTreeBuilder:
    def __init__(self, config, taxonomy):
        self.tree_nodes = {}
        self.leaf_count = {}
        self.taxonomy = taxonomy
        self.align_utils = AlignmentUtils(config)

    def prune_unifu_nodes(self, tree):
        for node in tree.traverse("preorder"):
            if len(node.children) == 1:
                node.delete()

    def add_tree_node(self, tree, nodeId, ranks, rank_level):
        if rank_level > 0:
            parent_level = rank_level            
            while ranks[parent_level] == "":
                parent_level -= 1
            parentId = ranks[parent_level]
        else:
            parentId = "root"

        if (parentId in self.tree_nodes):
            parentNode = self.tree_nodes[parentId]
        else:
            parentNode = self.add_tree_node(tree, parentId, ranks, parent_level-1)
            self.tree_nodes[parentId] = parentNode;

        newNode = parentNode.add_child()
        newNode.add_feature("name", nodeId)
        return newNode        

    def build(self, min_rank=0, max_seqs_per_leaf=1e9, clades_to_include=[], clades_to_ignore=[]):

        print "Number of nodes: %d" % self.taxonomy.seq_count()
        t0 = Tree()
        t0.add_feature("name", "root")
        self.tree_nodes["root"] = t0;
        self.leaf_count["root"] = 0;
        k = 0
        added = 0
        seq_ids = []
        # sequences are leafs of the tree, so they always have the lowest taxonomy level (e.g. "species"+1)        
        tax_seq_level = self.taxonomy.max_rank_level() + 1
        for sid, ranks in self.taxonomy.items():
            k += 1
            if k % 1000 == 0:
                print "Processed nodes: ", k, ", added: ", added, ", skipped: ", k - added

            # filter by minimum rank level            
            if ranks[min_rank] == "":
                continue       
    
            # filter by rank contraints (e.g. class Clostridia only)
            clade_is_ok = False

            # check against the inclusion list            
            if len(clades_to_include) > 0:
                for (rank_level, rank_name) in clades_to_include:            
                    if ranks[rank_level] == rank_name:
                        clade_is_ok = True
                        break
            else: # default: include all
                clade_is_ok = True

            # if sequence is about to be included, check it against the ignore list
            if clade_is_ok:
                for (rank_level, rank_name) in clades_to_ignore:
                    if ranks[rank_level] == rank_name:
                        clade_is_ok = False
                        break

            # final decision
            if not clade_is_ok:
                continue

            parent_level = tax_seq_level - 1            
            while ranks[parent_level] == "":
                parent_level -= 1
            parent_name = ranks[parent_level]
            if parent_name in self.tree_nodes:
                parent_node = self.tree_nodes[parent_name]
                # filter by max number of seqs (threshold depends from rank level, 
                # i.e. for genus there can be more seqs than for species)
                max_seq_per_rank = max_seqs_per_leaf * (tax_seq_level - parent_level)                
                if parent_name in self.leaf_count and self.leaf_count[parent_name] >= max_seq_per_rank:
                    continue

                old_sid_list = []
                for node in parent_node.children:
                    if node.is_leaf():
                        old_sid_list += [int(node.name)]
            else:
                old_sid_list = []

            # filter non-unique and invalid (e.g. "unaligned") sequences
#            if not self.align_utils.is_unique_sequence(old_sid_list, int(sid)):
#                continue

            if parent_name in self.leaf_count:
                self.leaf_count[parent_name] += 1
            else:
                # it'll be the first seq for a node, so init counter with 1                
                self.leaf_count[parent_name] = 1

            # all checks succeeded: add the sequence to the tree
            self.add_tree_node(t0, sid, ranks, parent_level)
            seq_ids += [sid]
            added += 1

        print "Total nodes in resulting tree: ", added

        self.prune_unifu_nodes(t0)
        return t0, seq_ids

    def close(self):
        self.align_utils.close()

class Taxonomy:
    def __init__(self):
        tree_nodes = []

    @staticmethod    
    def lineage_str(ranks):
        return ";".join(ranks).strip(';')

    @staticmethod    
    def lowest_assigned_rank_level(ranks):
        rank_level = len(ranks)-1
        while ranks[rank_level] == "":
            rank_level -= 1
        return rank_level

    @staticmethod    
    def lowest_assigned_rank(ranks):
        rank_level = Taxonomy.lowest_assigned_rank_level(ranks)
        return ranks[rank_level]

    def get_seq_ranks(self, seq_id):
        return []

    def seq_count(self):
        return 0

    def max_rank_level(self):
        return 0

    def items(self):
        return 0

class GGTaxonomyFile(Taxonomy):
    rank_placeholders = ["k__", "p__", "c__", "o__", "f__", "g__", "s__"]    

    def __init__(self, tax_fname):
        self.tax_fname = tax_fname        
        self.seq_ranks_map = {}
        self.load_taxonomy()
    
    def get_all_rank_names_at_level(self, rank_level):
        names = set([])
        for rks in self.seq_ranks_map.values():
            names.add(rks[rank_level])
        return names
    
    def extract_one_rank(self, rank_level, rank_name):
        target = {}
        remain = {}
        for item in self.items():
            sid = item[0]
            rks = item[1]
            if rks[rank_level] == rank_name:
                target[sid] = rks
            else:
                remain[sid] = rks
        return target, remain
    
    def seq_count(self):
        return len(self.seq_ranks_map)

    # zero-based! 
    def max_rank_level(self):
        return 6   # level of species in standard 7-level taxonomy
    
    def rank_level_name(self, rank_level):
        return { 0: "Kingdom",
                 1: "Phylum",
                 2: "Class",
                 3: "Order",
                 4: "Family",
                 5: "Genus",
                 6: "Species"
                }[rank_level]

    def get_seq_ranks(self, seq_id):
        return self.seq_ranks_map[seq_id]

    def lineage_str(self, seq_id, use_placeholders=False):
        ranks = list(self.seq_ranks_map[seq_id])
        if use_placeholders:        
            for i in range(len(ranks)):
                if ranks[i] == "":
                    ranks[i] = GGTaxonomyFile.rank_placeholders[i]
        return Taxonomy.lineage_str(ranks)

    def lowest_assigned_rank_level(self, seq_id):
        ranks = self.seq_ranks_map[seq_id]
        return Taxonomy.lowest_assigned_rank_level(ranks)

    def items(self):
        return self.seq_ranks_map.items()

    def load_taxonomy(self):
        print "Loading the taxonomy file..."
        fin = open(self.tax_fname)
        for line in fin:
            line = line.strip()
            toks = line.split("\t")
            sid = toks[0]
            ranks_str = toks[1]
            ranks = ranks_str.split(";")
            for i in range(len(ranks)):
                rank_name = ranks[i].strip()
                if rank_name in GGTaxonomyFile.rank_placeholders:
                    rank_name = ""
                ranks[i] = rank_name
                
            if len(ranks) < 7:
                ranks += [""] * (7 - len(ranks))
                print "WARNING: sequence " + sid + " has incomplete taxonomic annotation. Missing ranks were padded with empty values."  
            self.seq_ranks_map[sid] = ranks     

        fin.close()    

