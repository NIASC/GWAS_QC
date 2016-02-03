#!/usr/bin/python3.5
'''
Created on 15.12.2015
Only works on Atlas!
@author: paakkone
'''

__version__= "1.11"

import os
import sys
import subprocess
import argparse
import math
import numpy
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
        
#     def getPlinkSampleFileRow(self):
#         return "%s\t%s" % (self.indiv_id, self.fam_id)

class FailedSamples(object):
    """
    Stores the sample ID's for failed cases
    """
    def __init__(self):
        self.failed = []
        self.fail_types = ["sex", "miss_and_het", "ibd", "mds"]
        
    def addSample(self, indiv_id, fam_id, fail_type):
        new_sample = Sample(indiv_id, fam_id, fail_type)
        self.failed.append(new_sample)
        if fail_type not in self.fail_types:
            self.fail_types.append(fail_type)
            
    def writeFailedSamplesFile(self, out_fn, fail_type = None):
        # if fail type was given, check that it is one of the possible values
        if fail_type != None and fail_type not in self.fail_types:
            print ("Unknown fail type %s" % fail_type)
            print ("Allowed fail types: %s" % self.fail_types)
            raise ValueError(fail_type)
        # input ok, write the file
        indiv_ids_found = []
        with open(out_fn, 'w') as out_f:
            for sample in self.failed:
                indiv_id, fam_id, ftype = sample.getIDs()
                if indiv_id not in indiv_ids_found:
                    # new id, check optional fail type filtering
                    if fail_type == None or fail_type == ftype:
                        indiv_ids_found.append(indiv_id)
                        out_f.write("%s\t%s\n" % (indiv_id, fam_id))
        out_f.close()
        
    def getFailCount(self, fail_type):
        if fail_type not in self.fail_types:
            print ("Unknown fail type %s" % fail_type)
            print ("Allowed fail types: %s" % self.fail_types)
            raise ValueError(fail_type)
        fail_no = 0
        for sample in self.failed:
            if sample.getFailType() == fail_type:
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

    def getFailedIDs(self, limit):
        if type(limit) != type(1.0):
            print ("Wrong type for limit: %s %s" % (type(limit), limit))
            raise ValueError(limit)
        keys = self.samples.keys()
        failed_samples = []
        for key in keys:
            if float(self.samples[key]) >= limit:
                # failed sample
                fid, iid = key.split(self.delim)
                failed_samples.append((fid, iid))
        return failed_samples
                
    def getFIDAndIID(self, i):
        """
        Gets the index number and returns the corresponding FID and IID
        """
        keys = list(self.samples.keys())
        #print ("Keys[i] %s, keys[i].split(): %s" % (keys[i], keys[i].split(self.delim)))
        return keys[i].split(self.delim)
        
        
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
        # Check that output dir exists, or create it if it does not
        if not os.path.exists(self.args.output_dir):
            os.makedirs(self.args.output_dir)
            
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
        
    def runPlinkCommand(self, switches, in_fn = None, out_fn = None):
        """
        builds and runs plink command
        """
        cmd = self.buildPlinkCommand(switches, in_fn, out_fn)
        if cmd == None:
            print ("Could not run plink, exiting program")
            sys.exit()
        #print ("Run command %s" % cmd)
        self.runCommand(cmd)
        return True
        
    def runCommand(self, cmd):
        if os.name == "posix":
        # run SNPTest only in linux, otherwise unit tests will always fail
            print ("Running plink: %s" % cmd)
            p = subprocess.Popen(cmd.split(), stdin = subprocess.PIPE, stdout = subprocess.PIPE, close_fds = True)
            p.wait()
        
    def filterBySex(self, in_fn):
        """
        Steps 4-6
        """
        fail_type = "sex"
        root, head = os.path.split(in_fn)
        check_fn = head + "_%s" % fail_type
        check_fn = os.path.join(self.args.output_dir, check_fn)
        self.runPlinkCommand(["--check-sex"], in_fn, check_fn)
        # find problematic individuals
        fn = "%s.sexcheck" % check_fn
        #print ("in_fn: %s, \ncheck_fn: %s\nfn:%s" % (in_fn, check_fn, fn))
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
        #print ("No failed samples: %d" % len(self.failed_samples))
        return self.writeRemoveFileAndRemoveIndividuals(fail_type, in_fn, use_original = True)

    
    def writeRemoveFileAndRemoveIndividuals(self, fail_type, output_fn, use_original = False):
        fn = "remove_%s.txt" % fail_type
        fail_fn = os.path.join(self.args.output_dir, fn)
        self.writeFailedSamplesFile(fail_fn, fail_type = fail_type)
        root, head = os.path.split(output_fn)
        if use_original:
            input_fn = output_fn
        else:
            input_fn = os.path.join(self.args.output_dir, head)
        clean_set_fn = os.path.join(self.args.output_dir, head)
        if clean_set_fn.find("removed") == -1:
            clean_set_fn += "_removed"
        clean_set_fn += "_%s" % fail_type
        self.removeIndividualsFailingQC(fail_fn, input_fn, clean_set_fn)
        return clean_set_fn
    
    def removeIndividualsFailingQC(self, fail_fn, input_fn, output_fn):
        switch1 = "--remove %s" % fail_fn
        self.runPlinkCommand([switch1, "--make-bed"], input_fn, output_fn)
        return output_fn

    def filterHeterozygosityAndMissingRate(self, in_fn):
        """
        Steps 7 - 10
        """
        fail_type = "miss_and_het"
        out_fn = in_fn + "_missing" 
        # calculate missingness
        self.runPlinkCommand(["--missing","--allow-no-sex"], in_fn, out_fn)
        # store the results
        fn = "%s.imiss" % (out_fn)
        self.missingness = IMissFile(fn, self.delim)
        fmiss_vals = self.missingness.getMissingnessValues()
        # calculate heterozygosity
        out_fn = in_fn + "_heterozygosity"
        self.runPlinkCommand(["--het", "--allow-no-sex"], in_fn, out_fn)
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
                        indexes_to_delete.append(i)
                    else:
                        obs_het = (non_miss_no - hom) / (non_miss_no * 1.0)
                        het_vals.append(obs_het)
                        het_ids.append((fam_id, p_id))
        if len(indexes_to_delete) > 0:
            print ("Removing %d fmiss values due to zero values in heterozygosity" % len(indexes_to_delete))
            for ind in sorted(indexes_to_delete, reverse=True):
                del fmiss_vals[ind]
        # remove failed missingness values
        fmiss_limit = float(self.args.missingness)
        failed_fmiss = self.missingness.getFailedIDs(self.args.missingness)
        for fam_id, p_id in failed_fmiss:
            self.failed_samples.addSample(fam_id, p_id, fail_type)
            
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
        # plot the figure
        figure_title = "Heterozygosity rate and missingness"
        x_title = "Proportion of missing genotypes"
        y_title = "Heterozygosity rate"
        plotter = Plotter()
        col_map = plotter.plotScatterPlot(out_fn, fmiss_vals, het_vals, x_val_limit = fmiss_limit, 
                                          y_std_limit = std, figure_title = figure_title, 
                                          x_title = x_title, y_title = y_title, x_log_scale = True)
        return self.writeRemoveFileAndRemoveIndividuals(fail_type, in_fn)
        
    def filterDuplicates(self, in_fn):
        """
        Steps 11 - 13
        """
        fail_type = "ibd"
        # step 11 - filter high LD regions
        out_fn = in_fn + "_pruning"
        high_ldl_fn = os.path.join(self.dataset_path, "high-LD-regions.txt")
        switch1 = "--exclude %s" % high_ldl_fn
        self.runPlinkCommand([switch1,"--range", "--indep-pairwise 50 5 0.2"], in_fn, out_fn)
        # step 12 - generate pairwise IBS
        extract_str = "--extract %s.prune.in" % out_fn
        out_fn = in_fn + "_genome"
        self.runPlinkCommand([extract_str, "--genome"], in_fn, out_fn)
        # read .imiss and .genome files and identify those who have pi_hat over given limit
        genome_fn = "%s.genome" % out_fn
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
        #print ("HEADER ROW: %s" % header_row)
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
        return clean_set_fn, genome_fn
        
    def writeFailedSamplesFile(self, fn, fail_type = None):
        self.failed_samples.writeFailedSamplesFile(fn, fail_type)
        
    def getSTDLimits(self, c1_vals, c2_vals, c1_std, c2_std):
        mean1 = numpy.mean(c1_vals)
        std1 = numpy.std(c1_vals, ddof=1)
        upper_c1 = mean1 + (c1_std * std1)
        
        mean2 = numpy.mean(c2_vals)
        std2 = numpy.std(c2_vals, ddof=1)
        lower_c2 = mean2 - (c2_std * std2)
        return upper_c1, lower_c2
    
    def filterMultiDimensionalScaling(self, input_fn, gen_fn):
        """
        Step 17
        """
        fail_type = "mds"
        switch1 = "--read-genome %s" % gen_fn
        output_fn = input_fn + "_%s" % fail_type
        self.runPlinkCommand([switch1, "--cluster", "--mds-plot 10 "], input_fn, output_fn)
        # read and plot first two components
        mds_fn = output_fn + ".mds"
        header = True
        c1_vals = []
        c2_vals = []
        c1_limit = float(self.args.C1)
        c2_limit = float(self.args.C2)
        
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
                    c1_vals.append(c1)
                    c2_vals.append(c2)
                    mds_data.append([fid, iid, c1, c2])
        
        c1_std_upper_limit, c2_std_lower_limit = self.getSTDLimits(c1_vals, c2_vals, c1_limit, c2_limit)            
                    
        for fid, iid, c1, c2 in mds_data:
            if c1_limit > 0:
                if c1 >= c1_std_upper_limit:
                    #print ("C1 >= limit:IID: %s %3.3f >= %3.3f" % (iid, c1, c1_limit))
                    self.failed_samples.addSample(fid, iid, fail_type)
            if c2_limit > 0:
                if c2 <= c2_std_lower_limit:
                    #print ("C2 <= limit:IID: %s %3.3f <= %3.3f" % (iid, c2, c2_limit))
                    #print (line)
                    self.failed_samples.addSample(fid, iid, fail_type)
        # values gathered, make scatter plot
        figure_title = "Principal components C1 and C2"
        x_title = "C1"
        y_title = "C2"
        #std = 0.1
        out_fn = os.path.join(self.args.output_dir, "C1_and_C2_PCA_values.png")
        if c1_limit <= 1:
            c1_limit = None
        if c2_limit <= 1:
            c2_limit = None
        plotter = Plotter()
        col_map = plotter.plotScatterPlot(out_fn, c1_vals, c2_vals, x_std_limit = c1_limit, y_std_limit = c2_limit, 
                                          figure_title = figure_title, x_title = x_title,
                                          y_title = y_title, plot_std_y_axis = True)
#         col_map = plotter.plotScatterPlot(out_fn, c1_vals, c2_vals, y_val_limit= c2_limit, 
#                                           figure_title = figure_title, x_title = x_title,
#                                           y_title = y_title, plot_std_y_axis = True)
        return self.writeRemoveFileAndRemoveIndividuals(fail_type, input_fn)
    
#     def filterMissingDataRate(self, input_fn):
#         """
#         Filtering for excessive missing data rate (step 22) for the cleaned data set
#         """
#         out_fn = input_fn + "_missing"
#         self.runPlinkCommand(["--missing"], input_fn, out_fn)
#         fmiss_vals = []
#         fn = "%s.lmiss" % out_fn
#         with open(fn) as in_f:
#             for i, line in enumerate(in_f):
#                 if i > 0:
#                     # skipping header line
#                     words = line.split()
#                     fmiss_vals.append(float(words[4]))
#         self.plotter.plotHistogram(fmiss_vals)
        
    def filterSNPExclusion(self, input_fn):
        out_fn = input_fn + "_QCed"
        switch1 = "--geno %s" % self.args.geno
        switch2 = "--hwe %s" % self.args.hwe
        self.runPlinkCommand([switch1, switch2, "--make-bed"], input_fn, out_fn)
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
            no_removed = self.failed_samples.getFailCount("mds")
            out_f.write("Removed %d samples in multidimensional scaling check\n" % (no_removed))
            out_f.write("Call rate limit:\t%3.5f\n" % self.args.geno)
            out_f.write("HWE limit:\t%3.5f\n" % self.args.hwe)
            out_f.write("MAC limit for imputation:\t%3.5f\n" % self.args.mac)
    
    def main(self):
        cleaned_fn = self.filterBySex(self.base_name)
        cleaned_fn = self.filterHeterozygosityAndMissingRate(cleaned_fn)
        cleaned_fn, genome_fn = self.filterDuplicates(cleaned_fn)
        cleaned_fn = self.filterMultiDimensionalScaling(cleaned_fn, genome_fn)
        #cleaned_fn = self.filterMissingDataRate(cleaned_fn)
        cleaned_fn = self.filterSNPExclusion(cleaned_fn)
        cleaned_fn = self.filterForImputation(cleaned_fn)
        self.writeAnalysisReport()
        
def parseArguments(args):
    """
    Handles the argument parsing
    """
    parser = argparse.ArgumentParser(description = "Run QC filtering on GWAS datasets")
    parser.add_argument("--input", help = "Path to the dataset to be filtered")
    parser.add_argument("--output_dir", help = "Directory for the output files", required = True)
    parser.add_argument("--missingness", help = "Missingness limit, default is 0.05", type = float, default = 0.05)
    parser.add_argument("--heterozygosity", help = "Heterozygosity limit, default is 4.0 * sd", type = float, default = 4)
    parser.add_argument("--C1", help = "Principal component 1, default is sd * 4", type = float, default = 4)
    parser.add_argument("--C2", help = "Principal component 2, default is sd * 4", type = float, default = 4)
    #parser.add_argument("--C2", help = "Principal component 2, default is -0.06", type = float, default = -0.06)
    parser.add_argument("--pi_hat", help = "Pi-hat value of the genome file, default is 0.1", type = float, default = 0.1)
    parser.add_argument("--hwe", help = "HWE limit, default is 0.000001", type = float, default = 0.000001)
    parser.add_argument("--geno", help = "Call rate limit, default is 0.02", type = float, default = 0.02)
    parser.add_argument("--mac", help = "MAC treshold for imputation, default is 3", type = float, default = 3)
    parser.add_argument("--repeat", help = "Removing duplicates who's entry is repeated more than this limit (due to contamination). Default is 15",
                        type = int, default = 15)
    parser.add_argument("--plink_path", help = "Path to the plink executable", default = "/apps/genetics/bin/plink-1.9")
    parser.add_argument("--report_fn", help = "Analysis report file name", default = "analysis_report.txt")
    return parser.parse_args(args)

if __name__=='__main__':
    app = GWAS_QC()
    app.main()
    