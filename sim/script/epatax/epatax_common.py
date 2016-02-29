#!/usr/bin/env python

import os
import sys
import glob
import shutil
import datetime
from subprocess import call,STDOUT
import ConfigParser

class FileUtils:

    @staticmethod    
    def normalize_dir(dir_str):
        if dir_str and not dir_str.endswith("/"):
           dir_str += "/"
        return dir_str

    @staticmethod    
    def remove_if_exists(fname):
        try:
            os.remove(fname)
        except:
            pass

    @staticmethod    
    def rebase(fname, old_basedir, new_basedir):
        return fname.replace(old_basedir, new_basedir)    

class EpataxConfig:
    F_ALIGN_REF=1
    F_ALIGN_LBLQ=2
    F_ALIGN_QUERY=3

    F_TREE_MULTIF=11
    F_TREE_MULTIF_TAX=12
    F_TREE_OUTGR=13
    F_TREE_BIF_UNROOTED=14
    F_TREE_BIF_UNROOTED_LBL=15
    F_TREE_BIF_ROOTED_LBL=16
    F_TREE_BIF_ROOTED=17
    F_TREE_BIF_ROOTED_TAX=18
    F_TREE_EPA_RESULT=19
    F_TREE_EPA_RESULT_TAX=20

    F_SEQLIST_REF=31
    F_SEQLIST_LBLQ=32

    F_TAXONOMY_REF=41

    F_MAP_BRANCH_RANK=51
    F_OPT_MODEL=52

    F_QSUB_SCRIPT = 61,

    F_RESULT_TAX_ASSIGN=71
    F_RESULT_STATS=72
    F_RESULT_MIS=73
    F_RESULT_BRANCH_ASSIGN=74
    F_RESULT_LH_WEIGHTS=75

    def __init__(self, config_fname, results_name=""): 
        self.epatax_home = FileUtils.normalize_dir(os.path.abspath(""))
        self.raxml_outdir = "raxml_output/"
        self.raxml_outdir_abs = os.path.abspath(self.raxml_outdir)
        self.data_dir = FileUtils.normalize_dir(os.path.abspath("testdata/"))
        self.reftree_home = FileUtils.normalize_dir(os.path.abspath("reftree/"))
        self.temp_dir = FileUtils.normalize_dir(os.path.abspath("temp/"))
        self.results_home = FileUtils.normalize_dir(os.path.abspath("results/"))
        self.read_from_file(config_fname)
        if not results_name:
            results_name = self.reftree_name + "_" + datetime.datetime.now().strftime("%Y%m%d_%H%M")
        self.results_dir = self.results_home + results_name + "/"

    def read_from_file(self, config_fname):
        if not os.path.exists(config_fname):
            print "Config file not found: " + config_fname
            sys.exit()

        parser = ConfigParser.ConfigParser()
        parser.read(config_fname)
        self.raxml_remote_host = parser.get("raxml", "raxml_remote_host", "")
        self.raxml_home = FileUtils.normalize_dir(parser.get("raxml", "raxml_home"))
        self.raxml_exec = parser.get("raxml", "raxml_exec")
        self.raxml_exec_full = self.raxml_home + self.raxml_exec
        if self.raxml_remote_host in ["", "localhost"]:
            self.raxml_remote_call = False
            if not os.path.isdir(self.raxml_home):
                print "RAxML home directory not found: %s" % self.raxml_home
                sys.exit()
            if not os.path.isfile(self.raxml_exec_full):
                print "RAxML executable not found: %s" % self.raxml_exec_full
                sys.exit()
        else:
            self.raxml_remote_call = True

        self.raxml_model = parser.get("raxml", "raxml_model")
        self.raxml_threads = parser.get("raxml", "raxml_threads")
        self.raxml_cmd = [self.raxml_exec_full, "-p", "12345", "-T", self.raxml_threads, "-w", self.raxml_outdir_abs]

        self.epa_use_heuristic = parser.getboolean("raxml", "epa_use_heuristic")
        self.epa_heur_rate = parser.getfloat("raxml", "epa_heur_rate")
        self.epa_load_optmod = parser.getboolean("raxml", "epa_load_optmod")

        try:        
            self.run_on_cluster = parser.getboolean("cluster", "run_on_cluster")
            self.cluster_epatax_home = FileUtils.normalize_dir(parser.get("cluster", "cluster_epatax_home"))
            self.cluster_qsub_script = parser.get("cluster", "cluster_qsub_script")
        except ConfigParser.NoOptionError:
            self.run_on_cluster = False

        self.reftree_name = parser.get("reftree", "reftree_name")
        self.reftree_dir = self.reftree_home + self.reftree_name + "/"
        self.taxonomy_fname = parser.get("reftree", "taxonomy_file")

        self.reftree_min_rank = parser.getint("reftree", "min_rank")
        self.reftree_max_seqs_per_leaf = parser.getint("reftree", "max_seqs_per_leaf")
        clades_str = parser.get("reftree", "clades_to_include")
        self.reftree_clades_to_include = self.parse_clades(clades_str)
        try:
            clades_str = parser.get("reftree", "clades_to_ignore", "")
        except ConfigParser.NoOptionError:
            clades_str = ""            
        self.reftree_clades_to_ignore = self.parse_clades(clades_str)
            
        self.eval_taxonomy_fname = parser.get("eval", "taxonomy_file")
        try:
            clades_str = parser.get("eval", "ignored_clades", "")
        except ConfigParser.NoOptionError:
            clades_str = ""            
        self.eval_ignored_clades = self.parse_clades(clades_str)

        self.min_confidence = parser.getfloat("assignment", "min_confidence")

        self.db_host = parser.get("mysql", "db_host")
        self.db_user = parser.get("mysql", "db_user")
        self.db_pwd = parser.get("mysql", "db_pwd")
        self.db_name = parser.get("mysql", "db_name")

    def parse_clades(self, clades_str):
        clade_list = []
        try:        
            if clades_str:
                clades = clades_str.split(",")
                for clade in clades:
                    toks = clade.split("|")
                    clade_list += [(int(toks[0]), toks[1])]
        except:
            print "Invalid format in config parameter: clades_to_include"
            sys.exit()

        return clade_list

    def get_fname(self, file_type):
        return {
            EpataxConfig.F_ALIGN_REF: self.reftree_dir + self.reftree_name + ".fasta",
            EpataxConfig.F_ALIGN_LBLQ: self.temp_dir + self.reftree_name + "_lblq.fasta",
            EpataxConfig.F_ALIGN_QUERY: self.temp_dir + self.reftree_name + "_query.fasta",

            EpataxConfig.F_TREE_MULTIF: self.temp_dir + self.reftree_name + "_mfu.tre",
            EpataxConfig.F_TREE_MULTIF_TAX: self.temp_dir + self.reftree_name + "_mfu_tax.tre",
            EpataxConfig.F_TREE_OUTGR: self.temp_dir + self.reftree_name + "_outgr.tre",
            EpataxConfig.F_TREE_BIF_UNROOTED: self.reftree_dir + self.reftree_name + ".tre",
            EpataxConfig.F_TREE_BIF_UNROOTED_LBL: self.temp_dir + self.reftree_name + "_bfu_lbl.tre",
            EpataxConfig.F_TREE_BIF_ROOTED_LBL: self.temp_dir + self.reftree_name + "_bfr_lbl.tre",
            EpataxConfig.F_TREE_BIF_ROOTED: self.temp_dir + self.reftree_name + "_bfr.tre",
            EpataxConfig.F_TREE_BIF_ROOTED_TAX: self.temp_dir + self.reftree_name + "_bfr_tax.tre",
            EpataxConfig.F_TREE_EPA_RESULT: self.results_dir + self.reftree_name + "_epa.tre",
            EpataxConfig.F_TREE_EPA_RESULT_TAX: self.results_dir + self.reftree_name + "_epa_tax.tre",

            EpataxConfig.F_SEQLIST_REF: self.reftree_dir + self.reftree_name + "_seq.txt",
            EpataxConfig.F_SEQLIST_LBLQ: self.data_dir + "greengenes_lblq_seq.txt",

            EpataxConfig.F_TAXONOMY_REF: self.reftree_dir + self.reftree_name + "_tax.txt",

            EpataxConfig.F_MAP_BRANCH_RANK: self.reftree_dir + self.reftree_name + ".map",
            EpataxConfig.F_OPT_MODEL: self.reftree_dir + self.reftree_name + ".opt",

            EpataxConfig.F_QSUB_SCRIPT: self.temp_dir + self.reftree_name + "_sub.sh",

            EpataxConfig.F_RESULT_TAX_ASSIGN: self.results_dir + self.reftree_name + ".assigned.txt",
            EpataxConfig.F_RESULT_STATS: self.results_dir + self.reftree_name + ".stats.txt",
            EpataxConfig.F_RESULT_MIS: self.results_dir + self.reftree_name + ".mis.txt"
            }[file_type]

class RaxmlWrapper:

    def __init__(self, config): 
        self.config = config
    
    def make_raxml_fname(self, stem, job_name, absolute=True):
        fname = "RAxML_" + stem + "." + job_name
        if absolute:
            return self.config.raxml_outdir + fname
        else:
            return fname            

    def make_raxml_wildcard(self, job_name):
        return self.make_raxml_fname("*", job_name)

    def cleanup(self, job_name):
        raxml_out_mask = self.make_raxml_wildcard(job_name)
        for fl in glob.glob(raxml_out_mask):
            os.remove(fl)

    def reduce_alignment(self, align_fname, job_name="reduce"):
	    reduced_fname = align_fname + ".reduced"
	    FileUtils.remove_if_exists(reduced_fname)
	    self.run(job_name, ["-f", "c", "-s", align_fname])
	    self.cleanup(job_name)
	    return reduced_fname

    def run_epa(self, job_name, align_fname, reftree_fname, optmod_fname, silent=True):
        raxml_params = ["-f", "v", "-s", align_fname, "-t", reftree_fname]        
        if self.config.epa_use_heuristic:
            raxml_params += ["-G", str(self.config.epa_heur_rate)]
        if self.config.epa_load_optmod and optmod_fname:
            if os.path.isfile(optmod_fname):
                raxml_params += ["-R", optmod_fname]
            else:
                print "WARNING: Binary model file not found: %s" % optmod_fname
                print "WARNING: Model parameters will be estimated by RAxML"
                
        self.run(job_name, raxml_params, silent)

    def run(self, job_name, params, silent=True):
        self.cleanup(job_name)        
        params += ["-m", self.config.raxml_model, "-n", job_name]

        if self.config.run_on_cluster:
            self.run_cluster(params)
            return;        

        if self.config.raxml_remote_call:
            call_str = ["ssh", self.config.raxml_remote_host]
        else:
            call_str = []
        call_str += self.config.raxml_cmd + params
        if silent:        
            print ' '.join(call_str) + "\n"
            out_fname = self.make_raxml_fname("output", job_name)
            with open(out_fname, "w") as fout:
                call(call_str, stdout=fout, stderr=STDOUT)
        else:        
            call(call_str)

    def run_cluster(self, params):
        if self.config.raxml_remote_call:
            qsub_call_str = ["ssh", self.config.raxml_remote_host]
        else:
            qsub_call_str = []
        
        raxml_call_cmd = self.config.raxml_cmd + params        
        for i in range(len(raxml_call_cmd)):
            if isinstance(raxml_call_cmd[i], basestring):
                raxml_call_cmd[i] = FileUtils.rebase(raxml_call_cmd[i], self.config.epatax_home, self.config.cluster_epatax_home)
        raxml_call_str = ' '.join(raxml_call_cmd)
                
        script_fname = self.config.get_fname(EpataxConfig.F_QSUB_SCRIPT)
        FileUtils.remove_if_exists(script_fname)
        shutil.copy(self.config.cluster_qsub_script, script_fname)
        qsub_job_name = "epa"        
        with open(script_fname, "a") as fout:
            fout.write("#$ -N %s\n" % qsub_job_name)
            fout.write("\n")            
            fout.write(raxml_call_str + "\n")

        script_fname = FileUtils.rebase(script_fname, self.config.epatax_home, self.config.cluster_epatax_home)
        qsub_call_str += ["qsub", "-sync", "y", script_fname]

        print raxml_call_str + "\n"
        print ' '.join(qsub_call_str) + "\n"
#        sys.exit()

        call(qsub_call_str)

    def result_exists(self, job_name):
        if os.path.isfile(self.make_raxml_fname("result", job_name)):
            return True
        else:
            return False

    def epa_result_exists(self, job_name):
        if os.path.isfile(self.make_raxml_fname("labelledTree", job_name)):
            return True
        else:
            return False

    def copy_result_tree(self, job_name, dst_fname):
        src_fname = self.make_raxml_fname("result", job_name)
        shutil.copy(src_fname, dst_fname)

    def copy_optmod_params(self, job_name, dst_fname):
        src_fname = self.make_raxml_fname("binaryModelParameters", job_name)
        shutil.copy(src_fname, dst_fname)

    def copy_epa_orig_tree(self, job_name, dst_fname):
        src_fname = self.make_raxml_fname("originalLabelledTree", job_name)
        shutil.copy(src_fname, dst_fname)
        
    def copy_epa_result_tree(self, job_name, dst_fname):
        src_fname = self.make_raxml_fname("labelledTree", job_name)
        shutil.copy(src_fname, dst_fname)

