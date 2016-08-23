#!/apps/python3.5/anaconda3/bin/python
'''
#!/usr/bin/python3.5
Created on 15.12.2015
Only works on Atlas!
@author: paakkone
'''

__version__= "1.2"

import os
import sys
import subprocess
import argparse
import math
import numpy
import time
import random
from collections import OrderedDict
from plotter import Plotter


class Sample(object):
    """
    Single sample data
    """
    def __init__(self, indiv_id, fam_id, fail_type):
        self.indiv_id = indiv_id
        self.fam_id = fam_id
        self.fail_type = fail_type
    
    def getFailType(self):
        return self.fail_type
    
    def getIDs(self):
        return self.indiv_id, self.fam_id, self.fail_type
        
class FailedSamples(object):
    """
    Stores the sample ID's for failed cases
    """
    def __init__(self):
        self.failed = []
        self.fail_types = ["sex", "miss_and_het", "ibd", "mds"]
        
    def addSample(self, indiv_id, fam_id, fail_type):
        print ("Add sample: indiv_id %s, fam_id %s, fail_type %s" % (indiv_id, fam_id, fail_type))
        new_sample = Sample(indiv_id, fam_id, fail_type)
        self.failed.append(new_sample)
        if fail_type not in self.fail_types:
            self.fail_types.append(fail_type)
            
    def writeFailedSamplesFile(self, out_fn, fail_type = None):
        # if fail type was given, check that it is one of the possible values
#         if fail_type != None and fail_type not in self.fail_types:
#             print ("Unknown fail type %s" % fail_type)
#             print ("Allowed fail types: %s" % self.fail_types)
#             raise ValueError(fail_type)
        # input ok, write the file
        #print ("Writing failed samples file: %s  fail_type: %s" % (out_fn, fail_type))
        #print ("len(self.failed): %d" % len(self.failed))
        indiv_ids_found = []
        with open(out_fn, 'w') as out_f:
            #print ("Opened file %s" % out_fn)
            for sample in self.failed:
                indiv_id, fam_id, ftype = sample.getIDs()
                if indiv_id not in indiv_ids_found:
                    # new id, check optional fail type filtering
                    if fail_type == None or fail_type == ftype:
                        indiv_ids_found.append(indiv_id)
                        out_f.write("%s\t%s\n" % (indiv_id, fam_id))
                    # special case for mds, as previous round MDS extractions start with mdsrun1, mdsrun2
                    elif fail_type == "mdsruns":
                        if ftype.startswith("mdsrun"):
                            indiv_ids_found.append(indiv_id)
                            out_f.write("%s\t%s\n" % (indiv_id, fam_id))
                    #elif fail_type == "idb"
        out_f.close()
        #print ("Closed file %s" % out_fn)
        exists = os.path.exists(out_fn)
        #print ("Closed file exists: %s" % exists)
        
    def getFailCount(self, fail_type):
        if fail_type not in self.fail_types:
            print ("Unknown fail type %s" % fail_type)
            print ("Allowed fail types: %s" % self.fail_types)
            raise ValueError(fail_type)
        fail_no = 0
        for sample in self.failed:
            if sample.getFailType() == fail_type:
                fail_no += 1
            elif fail_type == "mds":
                sample_type = sample.getFailType()
                if sample_type.startswith("mds"):
                    fail_no += 1
        return fail_no
        
    def __getitem__(self, i):
        return self.failed[i]

    def __len__(self):
        return len(self.failed)

class IMissFile(object):
    """
    Class to read and hold imiss file data as that is needed in several places of the analysis
    """
    def __init__(self, fn, delim):
        self.fn = fn
        # stores family id, individual id and missingness value
        self.samples = OrderedDict()
        self.delim = delim
        self.readFile()
        
    def readFile(self):
        with open(self.fn) as in_f:
            for line in in_f:
                words = line.split()
                if words[0] != "FID":
                    key = words[0] + self.delim + words[1]
                    self.samples[key] = words[5]
    
    def getMissingnessValues(self, use_log = False):
        fmiss_vals = []
        for key in self.samples.keys():
            val = float(self.samples[key])
            if use_log:
                val = math.log10(val)
            fmiss_vals.append(val)
        return fmiss_vals
    
    def getMissingness(self, fid, iid):
        key = fid + self.delim + iid
        if key in self.samples.keys():
            return float(self.samples[key])
        # the fid iid pair not found from imiss file - raise error
        print ("No family id %s and individual id %s in file %s" % (fid, iid, self.fn))
        raise ValueError(key)

    def getFailedIDs(self, limit, get_indexes = False):
        if type(limit) != type(1.0):
            print ("Wrong type for limit: %s %s" % (type(limit), limit))
            raise ValueError(limit)
        print ("******************************")
        print ("Get failed IDs: limit: %s" % limit)
        keys = self.samples.keys()
        failed_samples = []
        for i, key in enumerate(keys):
            if float(self.samples[key]) >= limit:
                # failed sample
                if get_indexes:
                    failed_samples.append(i)
                    print ("failed value: %s, index: %d" % (float(self.samples[key]), i))
                else:
                    fid, iid = key.split(self.delim)
                    failed_samples.append((fid, iid))
                    print ("failed value: %s, index: %d  fid %s  iid %s" % (float(self.samples[key]), i, fid, iid))
        return failed_samples
                
    def writeFailedSamples(self, limit, out_fn):
        keys_to_del = []
        out_f = open(out_fn, 'w')
        for key in self.samples.keys():
            if float(self.samples[key]) >= limit:
                fid, iid = key.split(self.delim) 
                out_f.write("%s\t%s\n" % (fid, iid))
                keys_to_del.append(key)
        out_f.close()
        for key in keys_to_del:
            del self.samples[key]
            
    def getFIDAndIID(self, i):
        """
        Gets the index number and returns the corresponding FID and IID
        """
        keys = list(self.samples.keys())
        return keys[i].split(self.delim)
        
    def removeIndexes(self, indexes_to_delete):
        keys_to_delete = []
        for i, key in enumerate(self.samples.keys()):
            if i in indexes_to_delete:
                keys_to_delete.append(key)
        for key in keys_to_delete:
            del self.samples[key]
        print ("Removed %d indexes. Indexes: %s Keys: %s" % (len(indexes_to_delete), 
                                                             indexes_to_delete, keys_to_delete))
                
                
class GWAS_QC(object):
    '''
    classdocs
    '''
    def __init__(self, args = sys.argv[1:]):
        '''
        Constructor
        '''
        self.args = parseArguments(args)
        # get base path and dataset name
        root, head = os.path.split(self.args.input)
        self.dataset_path = root
        dot = head.find(".")
        if dot != -1:
            head = head[:dot]
        self.dataset_name = head
        self.base_name = os.path.join(self.dataset_path, self.dataset_name)
        # letter/string used to connect fid and iid to single key value
        self.delim = "*"
        # Failed samples
        self.failed_samples = FailedSamples()
        self.missingness = None
        self.global_geno = 0.02
        self.global_hwe = 0.000001
        # Check that output dir exists, or create it if it does not
        if not os.path.exists(self.args.output_dir):
            os.makedirs(self.args.output_dir)
        plink_runs_fn = os.path.join(self.args.output_dir, "plink_run_commands.log")
        self.plink_runs_f = open(plink_runs_fn, "w")
        
    def runPlinkCommand(self, switches, in_fn = None, out_fn = None):
        """
        builds and runs plink command
        """
        cmd = self.buildPlinkCommand(switches, in_fn, out_fn)
        if cmd == None:
            print ("Could not run plink, exiting program")
            sys.exit()
        self.runCommand(cmd)
        return True    
            
    def buildPlinkCommand(self, switches, in_fn = None, out_fn = None):
        """
        Gets switches for plink command and returns the command string
        """
        if out_fn == None:
            out_fn = os.path.join(self.args.output_dir, self.dataset_name)
        if in_fn == None:
            in_fn = self.args.input
        for ending in [".bim", ".bed", ".fam"]:
            fn = in_fn + ending
            if not os.path.exists(fn):
                print ("Plink input file %s does not exist!" % fn)
                return None
            if not os.path.isfile(fn):
                print ("Plink input file %s is not a file!" % fn)
                return None
        # start command and input file
        cmd = "%s --bfile %s" % (self.args.plink_path, in_fn)
        for sw in switches:
            cmd += " %s" % sw
        # output file
        cmd += " --out %s" % out_fn
        return cmd
        
    def runCommand(self, cmd):
        if os.name == "posix":
        # run SNPTest only in linux, otherwise unit tests will always fail
            print ("Running plink: %s" % cmd)
            self.plink_runs_f.write("%s\n" % cmd)
            p = subprocess.Popen(cmd.split(), stdin = subprocess.PIPE, stdout = subprocess.PIPE, close_fds = True)
            p.wait()
   
#     def writeFailedSamplesFile(self, fn, fail_type = None):
#         self.failed_samples.writeFailedSamplesFile(fn, fail_type)
     
    def writeRemoveFileAndRemoveIndividuals(self, fail_type, output_fn, use_original = False,
                                            replace_ending = 0):
        """
        Writes file which contains the samples which are removed, and then uses that file to 
        actually remove those samples. Returns the name of the resulting sample file.
        Normal behaviour is to extend the current file name with given file_type ending, but
        can also just replace last part (after last dash) with fail_type string.
        """
        # write the remove samples file
        fn = "remove_%s.txt" % fail_type
        print ("writeRemoveFileAndRemoveIndividuals: fail_type %s output_fn %s" % (fail_type, output_fn))
        fail_fn = os.path.join(self.args.output_dir, fn)
        self.failed_samples.writeFailedSamplesFile(fail_fn, fail_type)
        #self.writeFailedSamplesFile(fail_fn, fail_type = fail_type)
        root, head = os.path.split(output_fn)
        if use_original:
            input_fn = output_fn
        else:
            input_fn = os.path.join(self.args.output_dir, head)
        # make the name of the output file
        clean_set_fn = os.path.join(self.args.output_dir, head)
        if clean_set_fn.find("removed") == -1:
            clean_set_fn += "_removed"
        if replace_ending > 0:
            for i in range(replace_ending):
                last_dash = clean_set_fn.rfind("_")
                clean_set_fn = clean_set_fn[:last_dash]    
        clean_set_fn += "_%s" % fail_type
        # remove samples from the dataset
        print ("clean_set_fn: %s" % clean_set_fn)
        self.removeIndividualsFailingQC(fail_fn, input_fn, clean_set_fn)
        return clean_set_fn
    
    def removeIndividualsFailingQC(self, fail_fn, input_fn, output_fn):
        switch1 = "--remove %s" % fail_fn
        self.runPlinkCommand([switch1, "--make-bed"], input_fn, output_fn)
        return output_fn

    def filterBySex(self, in_fn):
        """
        Steps 4-6
        """
        fail_type = "sex"
        self.plink_runs_f.write("filterBySex, fail_type: %s\n" % fail_type)
        root, head = os.path.split(in_fn)
        check_fn = head + "_%s" % fail_type
        check_fn = os.path.join(self.args.output_dir, check_fn)
        self.runPlinkCommand(["--check-sex"], in_fn, check_fn)
        # find problematic individuals
        fn = "%s.sexcheck" % check_fn
        with open(fn) as in_f:
            for line in in_f:
                words = line.split()
                if words[4] == "PROBLEM":
                    sex1 = words[2]
                    sex2 = words[3]
                    if sex1 in ["1", "2"] and sex2 in ["1", "2"]:
                        fam_id = words[0]
                        p_id = words[1]
                        self.failed_samples.addSample(fam_id, p_id, fail_type)
        return self.writeRemoveFileAndRemoveIndividuals(fail_type, in_fn, use_original = True)

    def filterHeterozygosityAndMissingRate(self, in_fn):
        """
        Steps 7 - 10
        """
        fail_type = "miss_and_het"
        self.plink_runs_f.write("\nfilterHeterozygosityAndMissingRate, fail_type: %s\n" % fail_type)
        out_fn = in_fn + "_missing" 
        # calculate missingness
        self.runPlinkCommand(["--missing","--allow-no-sex"], in_fn, out_fn)
        # store the results
        fn = "%s.imiss" % (out_fn)
        self.missingness = IMissFile(fn, self.delim)
        #fmiss_vals = self.missingness.getMissingnessValues()
        #print ("len(miss_vals): %d" % (len(fmiss_vals)))
        # calculate heterozygosity
        out_fn = in_fn + "_heterozygosity"
        switch1 = "--geno %s" % self.global_geno
        #switch1 = "--geno 0.1"
        switch2 = "--hwe %s" % self.global_hwe
        #self.runPlinkCommand([switch1, switch2, "--het", "--allow-no-sex"], in_fn, out_fn)
        missing_remove_fn = os.path.join(self.args.output_dir, "remove_missingness_for_heterozygosity_test.txt")
        #self.failed_samples.writeFailedSamplesFile(missing_remove_fn, fail_type)
        self.missingness.writeFailedSamples(self.args.missingness, missing_remove_fn)
        fmiss_vals = self.missingness.getMissingnessValues()
        switch3 = "--remove %s" % missing_remove_fn
        self.runPlinkCommand([switch1, switch2, switch3, "--het", "--allow-no-sex"], in_fn, out_fn)
        fn = "%s.het" % (out_fn)
        het_vals = []
        het_ids = []
        indexes_to_delete = []
        with open(fn) as in_f:
            for i, line in enumerate(in_f):
                words = line.split()
                if words[0] != "FID":
                    fam_id = words[0]
                    p_id = words[1]
                    hom = float(words[2])
                    non_miss_no = float(words[4])
                    if non_miss_no == 0:
                        indexes_to_delete.append(i-1)
                        self.failed_samples.addSample(fam_id, p_id, fail_type)
                    else:
                        obs_het = (non_miss_no - hom) / (non_miss_no * 1.0)
                        het_vals.append(obs_het)
                        het_ids.append((fam_id, p_id))
        # first remove missingness values which had non_miss_no = 0. These won't be plotted (as there is no 
        # corresponding heterozygosity value)
        if len(indexes_to_delete) > 0:
            print ("Removing %d fmiss values due to zero values in heterozygosity" % len(indexes_to_delete))
            for ind in sorted(indexes_to_delete, reverse=True):
                print ("ind: %s len(fmiss_vals): %s" % (ind, len(fmiss_vals)))
                #print ("Deleting fmiss_vals(%s), het_ids[ind] (%s) of ind %d" % (fmiss_vals[ind], het_ids[ind], ind))
                del fmiss_vals[ind]
            self.missingness.removeIndexes(indexes_to_delete)
        fmiss_vals2 = self.missingness.getMissingnessValues()
        print ("len(miss_vals): %d len(fmiss_vals2): %d  len(het_vals): %d" % (len(fmiss_vals), len(fmiss_vals2), len(het_vals)))
        cm = ["green"] * len(het_vals)
        # remove failed missingness values
        fmiss_limit = float(self.args.missingness)
        failed_fmiss = self.missingness.getFailedIDs(self.args.missingness)
        failed_indexes = self.missingness.getFailedIDs(self.args.missingness, get_indexes = True)
        for fam_id, p_id in failed_fmiss:
            self.failed_samples.addSample(fam_id, p_id, fail_type)
        print ("\n***********************************************")
        for i in failed_indexes:
            cm[i] = "red"
            print ("Failed missingness: het_val:%s  fmiss_val:%s" % (het_vals[i], fmiss_vals[i]))
        print ("Number of failed fmiss: %d,  %d" % (len(failed_fmiss), len(failed_indexes)))
        
        if len(failed_fmiss) != len(failed_indexes):
            print ("Number of errors different!")
            raise ValueError
        print ("\n***********************************************")
        out_fn = os.path.join(self.args.output_dir, "missingness_and_heterozygosity.png")
        std = float(self.args.heterozygosity)
        # remove failed heterozygosity values
        mean_heterozygosity = numpy.mean(het_vals)
        std_heterozygosity = numpy.std(het_vals, ddof=1)
        upper = mean_heterozygosity + (std * std_heterozygosity)
        lower = mean_heterozygosity - (std * std_heterozygosity)
        for i in range(len(het_vals)):
            if het_vals[i] > upper or het_vals[i] < lower:
                fam_id, p_id = het_ids[i]
                self.failed_samples.addSample(fam_id, p_id, fail_type)
                cm[i] = "red"
                print ("Failed heterozygosity: het_val: %s, fmiss_val: %s" % (het_vals[i], fmiss_vals[i]))
        # plot the figure
        figure_title = "Heterozygosity rate and missingness"
        x_title = "Proportion of missing genotypes"
        y_title = "Heterozygosity rate"
        plotter = Plotter()
        plotter.plotScatterPlot(out_fn, fmiss_vals, het_vals, x_val_limit = fmiss_limit, 
                                  y_std_limit = std, figure_title = figure_title, 
                                  x_title = x_title, y_title = y_title, x_log_scale = True, col_map = cm)
#         col_map = plotter.plotScatterPlot(out_fn, fmiss_vals, het_vals, x_val_limit = fmiss_limit, 
#                                           y_std_limit = std, figure_title = figure_title, 
#                                           x_title = x_title, y_title = y_title, x_log_scale = True, col_map = cm)
#         print ("\n***********************************************")
#         print ("Color map from routine: %d Reds: %d   Greens: %d" % (len(cm), cm.count("red"), cm.count("green")))
#         print ("Color map from plotter: %d Reds: %d   Greens: %d" % (len(col_map), col_map.count("red"), col_map.count("green")))
#         print ("\n***********************************************")
        return self.writeRemoveFileAndRemoveIndividuals(fail_type, in_fn)
        
    def filterDuplicates(self, in_fn, only_genome = False):
        """
        Steps 11 - 13
        """
        fail_type = "ibd"
        if only_genome:
            self.plink_runs_f.write("\n***MDSLoop** filterDuplicates, fail_type: %s\n" % fail_type)
        else:
            self.plink_runs_f.write("\nfilterDuplicates, fail_type: %s\n" % fail_type)
        # step 11 - filter high LD regions
        out_fn = in_fn + "_pruning"
        
        print ("********************************")
        print ("filterDuplicates: in_fn: %s" % in_fn)
        print ("only_genome: %s out_fn: %s" % (only_genome, out_fn))
        
        if not os.path.isfile(self.args.high_ld_regions):
            print ("Could not find file %s for duplicate filtering" % self.args.high_ld_regions)
            sys.exit()
        switch1 = "--geno %s" % self.global_geno
        switch2 = "--hwe %s" % self.global_hwe
        switch3 = "--exclude %s --range" % self.args.high_ld_regions
        self.runPlinkCommand([switch1,switch2, switch3,"--indep-pairwise 50 5 0.2"], in_fn, out_fn)
        #self.runPlinkCommand([switch3, "--indep-pairwise 50 5 0.2"], in_fn, out_fn)
        # step 12 - generate pairwise IBS
        extract_str = "--extract %s.prune.in" % out_fn
        out_fn = in_fn + "_genome"
        self.runPlinkCommand([extract_str, "--genome"], in_fn, out_fn)
        # read .imiss and .genome files and identify those who have pi_hat over given limit
        genome_fn = "%s.genome" % out_fn
        if only_genome:
            clean_set_fn = self.writeRemoveFileAndRemoveIndividuals(fail_type, in_fn)
            return clean_set_fn, genome_fn
        remove = []
        key_counter = {}
        repeat_limit = int(self.args.repeat)
        pihat_over_limit = []
        header_row = []
        with open(genome_fn) as in_f:
            for line in in_f:
                words = line.split()
                if len(header_row) == 0:
                    # skip header
                    header_row = [words[0],words[1],words[2],words[3],words[9]]
                else:
                    fid1 = words[0]
                    iid1 = words[1]
                    fid2 = words[2]
                    iid2 = words[3]
                    key1 = fid1 + self.delim + iid1
                    key2 = fid2 + self.delim + iid2
                    pi_hat = float(words[9])
                    if pi_hat > float(self.args.pi_hat):
                        pihat_over_limit.append((key1, key2, words[9]))
                        for key in [key1, key2]:
                            if key not in key_counter.keys():
                                key_counter[key] = 1
                            else:
                                key_counter[key] += 1
                            if key_counter[key] > repeat_limit:
                                if key not in remove:
                                    remove.append(key)
        # go through the pihat list again and remove those who will be removed
        fn = "PI_hat_over_%1.2f.txt" % self.args.pi_hat
        out_fn = os.path.join(self.args.output_dir, fn)
        h_row = ""
        with open(out_fn, 'w') as out_f:
            for header in header_row:
                h_row += "%s\t" % header
            out_f.write("%s\n" % h_row[:-1])
            for key1, key2, pi_hat in pihat_over_limit:
                if (key1 in remove) or (key2 in remove):
                    pass
                else:
                    fid1, pid1 = key1.split(self.delim)
                    fid2, pid2 = key2.split(self.delim)
                    out_f.write("%s\t%s\t%s\t%s\t%s\n" % (fid1, pid1,fid2, pid2, pi_hat))
        # remove contaminations (too many repeats)
        for ids in remove:
            fid, pid = ids.split(self.delim)
            self.failed_samples.addSample(fid, pid, fail_type)
        clean_set_fn = self.writeRemoveFileAndRemoveIndividuals(fail_type, in_fn)
        print ("********************************")
        return clean_set_fn, genome_fn, fn
        
        
    def getSTDLimits(self, c1_vals, c2_vals, c1_std, c2_std):
        mean1 = numpy.mean(c1_vals)
        std1 = numpy.std(c1_vals, ddof=1)
        upper_c1 = mean1 + (c1_std * std1)
        
        mean2 = numpy.mean(c2_vals)
        std2 = numpy.std(c2_vals, ddof=1)
        lower_c2 = mean2 - (c2_std * std2)
        upper_c2 = mean2 + (c2_std * std2)
        return upper_c1, lower_c2, upper_c2
    
    def addIDCount(self, id_counts, fid, iid):
        new_id = fid + "&" + iid
        if new_id not in id_counts.keys():
            id_counts[new_id] = 0
        id_counts[new_id] += 1
        return id_counts
    
    def getIDCounts(self, pihat_rows):
        """
        Counts the occurances of the ID's
        """
        print ("**********Get ID counts")
        print ("First element: %s" % pihat_rows[0])
        id_counts = {}
        for fid1, iid1, fid2, iid2, pi_hat in pihat_rows:
            id_counts = self.addIDCount(id_counts, fid1, iid1)
            id_counts = self.addIDCount(id_counts, fid2, iid2)
        return id_counts
    
    def removeRelatedFromMDS(self, pi_val_fn, dataset_fn):
        """
        Reads in the pihat value file which contains all sample pairs with pihat values greater than given limit.
        For MDS these samples need to be removed (or half of them). For that it checks for repeated ID's and
        removes those which are several times in the file first, until each sample is just once in the file.
        Then from those one sample per row is randomly selected for removal.
        """
        fail_type = "mds_rel"
        self.plink_runs_f.write("\nremoveRelatedFromMDS, fail_type: %s\n" % fail_type)
        pihat_rows = []
        # id list of id's which will be removed
        ids_to_remove = []
        pi_val_fn = os.path.join(self.args.output_dir, pi_val_fn)
        header = True
        with open(pi_val_fn) as in_f:
            for line in in_f:
                if header:
                    header = False
                else:
                    words = line.split()
                    pihat_rows.append(words)
                    if len(words) != 5:
                        print ("*********Only four entries! %s" % words)
        # values read, now go through rows and remove IDs from it starting with the most common id's
        done = False
        if len(pihat_rows) == 0:
            # skip this step if there were no pihat values over the treshold in the data
            done = True
        while not done:
            id_counts = self.getIDCounts(pihat_rows)
            max_repetitions = max(id_counts.values())
            if max_repetitions > 1:
                # find all max value ids
                max_val_sample_ids = [sid for sid in id_counts if id_counts[sid] == max_repetitions]
                # choose randomly one of the max repetition value ID's
                selected_id = random.choice(max_val_sample_ids)
            else:
                # go through rest of the lines and choose from each line randomly one key
                fid1, iid1, fid2, iid2, pi_hat = pihat_rows[0]
                selected_id = fid1 + "&" + iid1
                if random.randint(1,2) == 2:
                    selected_id = fid2 + "&" + iid2                    
            # put selected value to ids_to_remove list
            ids_to_remove.append(selected_id)
            # remove rows which contain selected id
            pihat_rows = self.removeSelectedIDRows(pihat_rows, selected_id)
            # check if there is any rows left
            if len(pihat_rows) == 0:
                done = True
        # write out to_remove list
        for sid in ids_to_remove:
            fid, iid = sid.split("&")
            self.failed_samples.addSample(fid, iid, fail_type)
        return self.writeRemoveFileAndRemoveIndividuals(fail_type, dataset_fn)
    
    def removeSelectedIDRows(self, pihat_rows, selected_id):
        """
        Delete rows from pihat_rows which contains the selected ID
        """
        remove_fid, remove_iid = selected_id.split("&")
        # don't actually delete items but create a new list with those rows which don't have the ID (simpler than deleting)
        remaining_rows = []
        #for i, fid1, iid1, fid2, iid2, pi_hat in enumerate(pihat_rows):
        for i, row in enumerate(pihat_rows):
            fid1, iid1, fid2, iid2, pi_hat = row
            if not (fid1 == remove_fid and iid1 == remove_iid) and not (fid2 == remove_fid and iid2 == remove_iid):
                remaining_rows.append(pihat_rows[i])
        return remaining_rows
    
    def getMSDLimit(self, round_no):
        string_limits = (self.args.C1_string, self.args.C2_string)
        float_limits = (self.args.C1, self.args.C2)
        limits = []
        for i in range(2):
            if string_limits[i] == "":
                limits.append(float(float_limits[i]))
            else:
                # limits string are _ separated list of floating point values
                fields = string_limits[i].split("_")
                field = fields[round_no - 1]
                limits.append(float(field))
        return limits
 
    def filterMultiDimensionalScaling(self, input_fn, gen_fn, round_no):
        """
        Step 17
        """
        fail_type = "mdsrun%d" % round_no
        self.plink_runs_f.write("\n***MDSLoop** filterMDS, fail_type: %s\n" % fail_type)
        #fail_type = "mdsruns"
        switch1 = "--read-genome %s" % gen_fn
        if round_no == 1:
            output_fn = input_fn + "_%s" % fail_type
        else:
            # remove _idb
            last_dash = input_fn.rfind("_")
            # remove previous rounds _mdsrun#
            last_dash = input_fn[:last_dash].rfind("_")
            output_fn = input_fn[:last_dash] + "_%s" % fail_type
            
        print ("********************************")
        print ("filterMDS: in_fn: %s" % input_fn)
        print ("round_no: %s   out_fn: %s" % (round_no, output_fn))
        
        self.runPlinkCommand([switch1, "--cluster", "--mds-plot 10 "], input_fn, output_fn)
        # read and plot first two components
        mds_fn = output_fn + ".mds"
        header = True
        c1_vals = []
        c2_vals = []
        c1_limit, c2_limit = self.getMSDLimit(round_no)
        #c1_limit = float(self.args.C1)
        #c2_limit = float(self.args.C2)
        
        mds_data = []
        with open(mds_fn) as mds_f:
            for line in mds_f:
                if header:
                    # skip header line
                    header= False
                else:
                    words = line.split()
                    fid = words[0]
                    iid = words[1]
                    c1 = float(words[3])
                    c2 = float(words[4])
                    if words[3] == "nan" or words[4] =="nan":
#                    if type(c1) != type(float) or type(c2) != type(float):
                        print ("*************************\n\n")
                        print ("Bad line: %s" % line)
                        print ("c1: %s   c2: %s" % (words[3], words[4]))
                        print ("\n*********")
                        sys.exit()
                    c1_vals.append(c1)
                    c2_vals.append(c2)
                    mds_data.append([fid, iid, c1, c2])
        
        c1_std_upper_limit, c2_std_lower_limit, c2_std_upper_limit = self.getSTDLimits(c1_vals, c2_vals, c1_limit, c2_limit)            
        cm = ["green"] * len(mds_data)     
        i = 0
        no_violations = 0
        for fid, iid, c1, c2 in mds_data:
            limit_violation = False
            if c1_limit > 0:
                if c1 >= c1_std_upper_limit:
                    limit_violation = True
            if c2_limit > 0:
                if c2 <= c2_std_lower_limit:
                    limit_violation = True
                if c2 >= c2_std_upper_limit:
                    limit_violation = True
            if limit_violation:
                self.failed_samples.addSample(fid, iid, fail_type)
                cm[i] = "red"
                no_violations += 1    
            i += 1
        # values gathered, make scatter plot
        figure_title = "Principal components C1 and C2"
        x_title = "C1"
        y_title = "C2"
        name = "C1_and_C2_PCA_values_round_%d.png" % round_no
        out_fn = os.path.join(self.args.output_dir, name)
        if c1_limit <= 1:
            c1_limit = None
        if c2_limit <= 1:
            c2_limit = None
        plotter = Plotter()
        print ("********************\n\n")
        print ("len(c1_vals) %d: len(c2_vals) %d" % (len(c1_vals), len(c2_vals)))
        print ("c1_limit: %s, c2_limit: %s" % (c1_limit, c2_limit))
        print ("len(col_map) %d" % len(cm))
        for i in range(len(c1_vals)):
            if type(c1_vals[i]) != type(float):
                print ("Non-float number in c1: %s" % c1_vals[i])
            if type(c2_vals[i]) != type(float):
                print ("Non-float number in c1: %s" % c2_vals[i])
                
        print ("\n\n**********************\n")
        plotter.plotScatterPlot(out_fn, c1_vals, c2_vals, x_std_limit = c1_limit, y_std_limit = c2_limit, 
                                  figure_title = figure_title, x_title = x_title,
                                  y_title = y_title, plot_std_y_axis = True, col_map = cm)
#         col_map = plotter.plotScatterPlot(out_fn, c1_vals, c2_vals, x_std_limit = c1_limit, y_std_limit = c2_limit, 
#                                           figure_title = figure_title, x_title = x_title,
#                                           y_title = y_title, plot_std_y_axis = True, col_map = cm)
        rep_no = 1
        if round_no > 1:
            rep_no = 2
        out_fn = self.writeRemoveFileAndRemoveIndividuals(fail_type, input_fn, replace_ending = rep_no)
        print ("No violations: %s" % no_violations)
        print ("********************************")
        return out_fn, no_violations
        #return self.writeRemoveFileAndRemoveIndividuals(fail_type, input_fn)
    
    def filterSNPExclusion(self, input_fn):
        self.plink_runs_f.write("\nfilterSNPExclusion\n")
        out_fn = input_fn + "_QCed"
        switch1 = "--geno %s" % self.args.geno
        switch2 = "--hwe %s" % self.args.hwe
        
        missing_remove_fn = os.path.join(self.args.output_dir, "remove_missingness_for_heterozygosity_test.txt")
        switch3 = "--remove %s" % missing_remove_fn
        self.runPlinkCommand([switch1, switch2, switch3, "--make-bed"], input_fn, out_fn)
        
        #self.runPlinkCommand([switch1, switch2, "--make-bed"], input_fn, out_fn)
        return out_fn
    
    def filterByGenotypeCallRates(self):
        """
        This step (24) is only for case-control cases. 
        """
        pass
    
    def filterForImputation(self, input_fn):
        """
        Use bim and freq file to filter non-autosomal snips and snps with 
        MAC value less than  treshold
        """
        freq_fn = input_fn + "_freq"
        self.plink_runs_f.write("\nfilterForImputation\n")
        self.runPlinkCommand(["--freq"], input_fn, freq_fn)
        filtered_variants = OrderedDict()
        # filter non autosomal or x variants (anything else but chromosomes 1-23)
        with open(input_fn + ".bim") as bim_f:
            for line in bim_f:
                words = line.split()
                remove = False
                try:
                    chrom = int(words[0])
                    if chrom < 1 or chrom > 23:
                        remove = True
                except ValueError:
                    # if words[0] was X or Y, int() will fail
                    remove = True
                if remove:
                    rsid = words[1]
                    if rsid not in filtered_variants.keys():
                        filtered_variants[rsid] = 0
        # filter low freq variants
        header = True
        with open(freq_fn + ".frq") as freq_f:
            for line in freq_f:
                if header:
                    header = False
                else:
                    words = line.split()
                    maf = float(words[4])
                    nchrobs = int(words[5])
                    mac = int(round(nchrobs * maf))
                    if mac < self.args.mac:
                        rsid = words[1]
                        if rsid not in filtered_variants.keys():
                            filtered_variants[rsid] = 0
        # write filtered rsid:s to file
        fail_fn = "Low_freq_and_Non_autosomal_or_X_SNPs.txt"
        fail_fn = os.path.join(self.args.output_dir, fail_fn)
        with open(fail_fn, 'w') as fail_f:
            for rsid in filtered_variants.keys():
                fail_f.write("%s\n" % (rsid))         
        # use plink to filter these snps
        out_fn = input_fn + "_for_imputation"
        switch1 = "--exclude %s" % fail_fn
        self.runPlinkCommand([switch1, "--make-bed"], input_fn, out_fn)
        return out_fn
        
    def writeAnalysisReport(self):
        report_fn = os.path.join(self.args.output_dir, self.args.report_fn)
        with open(report_fn, 'w') as out_f:
            out_f.write("Program version: %s\n" % __version__)
            out_f.write("Running date: %s\n\n" % time.ctime())
            out_f.write("Input dataset:\t%s\n" % self.args.input)
            no_removed = self.failed_samples.getFailCount("sex")
            out_f.write("Removed %d samples in sex check\n" % (no_removed))
            out_f.write("Missingness rate:\t%3.5f\n" % self.args.missingness)
            out_f.write("Heterozygosity rate:\t%3.5f * SD\n" % self.args.heterozygosity)
            no_removed = self.failed_samples.getFailCount("miss_and_het")
            out_f.write("Removed %d samples in missingness and heterotsygosity check\n" % (no_removed))
            out_f.write("Pi-hat value:\t%3.5f\n" % self.args.pi_hat)
            no_removed = self.failed_samples.getFailCount("ibd")
            out_f.write("Removed %d samples in contamination check\n" % (no_removed))
            out_f.write("Principal component 1 cutoff:\t%3.5f * SD\n" % self.args.C1)
            out_f.write("Principal component 2 cutoff:\t%3.5f * SD\n" % self.args.C2)
            out_f.write("Maximum number of MDS rounds: %d\n" % self.args.max_mds_rounds)
            out_f.write("Principal component 1 cutoff string: %s\n" % self.args.C1_string)
            out_f.write("Principal component 2 cutoff string: %s\n" % self.args.C2_string)
            no_removed = self.failed_samples.getFailCount("mds")
            out_f.write("Removed %d samples in multidimensional scaling check\n" % (no_removed))
            out_f.write("Call rate limit:\t%3.5f\n" % self.args.geno)
            out_f.write("HWE limit:\t%3.5f\n" % self.args.hwe)
            out_f.write("MAC limit for imputation:\t%3.5f\n" % self.args.mac)
    
    def main(self):
        cleaned_fn = self.filterBySex(self.base_name)
        cleaned_fn = self.filterHeterozygosityAndMissingRate(cleaned_fn)
        before_mds_fn, genome_fn, related_fn = self.filterDuplicates(cleaned_fn)
        
        final = False
        mds_fn = self.removeRelatedFromMDS(related_fn, before_mds_fn)
        round_no = 0
        while not final:
            round_no += 1
            # do MDS
            mds_fn, genome_fn = self.filterDuplicates(mds_fn, only_genome = True)
            mds_fn, no_overlimit = self.filterMultiDimensionalScaling(mds_fn, genome_fn, round_no)
            # redo genome for next round of MDS
            if (no_overlimit == 0) or (round_no > self.args.max_mds_rounds):
                final = True
                
        # MDS done, remove samples identified by MDS but keep related samples that were removed for MDS
        self.plink_runs_f.write("\nremove samples that failed in MDS runs\n")
        after_mds_fn = self.writeRemoveFileAndRemoveIndividuals("mdsruns", before_mds_fn)
        #after_mds_fn = self.removeSamplesByMDS(before_mds_fn)
            
        cleaned_fn = self.filterSNPExclusion(after_mds_fn)
        cleaned_fn = self.filterForImputation(cleaned_fn)
        self.writeAnalysisReport()
        
        self.plink_runs_f.close()
        
def parseArguments(args):
    """
    Handles the argument parsing
    """
    parser = argparse.ArgumentParser(prog = "gwas_qc",
                                     description = "Run QC filtering on GWAS datasets",
                                     formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--input", help = "Path to the dataset to be filtered")
    parser.add_argument("--output_dir", help = "Directory for the output files", required = True)
    parser.add_argument("--missingness", help = "Missingness limit", type = float, default = 0.05)
    parser.add_argument("--heterozygosity", help = "Heterozygosity limit = default * sd", type = float, default = 4)
    parser.add_argument("--C1", help = "Principal component 1 limit = default * sd", type = float, default = 4)
    parser.add_argument("--C2", help = "Principal component 2 limit = default * sd", type = float, default = 4)
    parser.add_argument("--C1_string", help = "Principal component 1 limit string = default * sd", default = "")
    parser.add_argument("--C2_string", help = "Principal component 2 limit string = default * sd", default = "")
    parser.add_argument("--max_mds_rounds", help = "Maximum number of MDS rounds", type = int, default = 5)
    parser.add_argument("--pi_hat", help = "Pi-hat limit for reporting samples", type = float, default = 0.1)
    parser.add_argument("--pi_hat_removal_limit", help = "Pi-hat limit for removing samples", type = float, default = 0.95)
    parser.add_argument("--hwe", help = "HWE limit", type = float, default = 0.000001)
    parser.add_argument("--geno", help = "Call rate limit", type = float, default = 0.02)
    parser.add_argument("--mac", help = "MAC treshold for imputation", type = float, default = 3)
    parser.add_argument("--repeat", help = "Removing duplicates who's entry is repeated more than this limit (due to contamination).",
                        type = int, default = 15)
    parser.add_argument("--plink_path", help = "Path to the plink executable", default = "/apps/genetics/bin/plink-1.9")
    parser.add_argument("--report_fn", help = "Analysis report file name", default = "analysis_report.txt")
    parser.add_argument("--high_ld_regions", help = "Path to the file used in duplicate filtering",
                        default = "/fs/projects/finngen/misc/high-LD-regions.txt")
    parser.add_argument("--version", help = "Shows the version number of the program",
                        action = 'version', version = "%(prog)s {version}".format(version=__version__))
    return parser.parse_args(args)

if __name__=='__main__':
    app = GWAS_QC()
    app.main()
    