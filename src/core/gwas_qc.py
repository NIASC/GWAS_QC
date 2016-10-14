#!/apps/python3.5/anaconda3/bin/python
'''
Created on 18.8.2016

@author: paakkone
'''
__version__= "1.2"

import sys
import os
import time
import argparse
from filterBySex import FilterBySex
from filterHetAndMissing import FilterHeterozygosityAndMissingRate
from IBSCalculation import IBSCalculation
from removeRelatedFromMDS import RemoveRelatedFromMDS
from filterMDS import FilterMDS
from filterSNPExclusion import FilterSNPExclusion
from filterForImputation import FilterForImputation

class GWAS_QC(object):
    def __init__(self, args = sys.argv[1:]):
        '''
        Constructor
        '''
        self.args = parseArguments(args)
        if not os.path.exists(self.args.output_dir):
            os.makedirs(self.args.output_dir)
        plink_fn = "plink_commands.txt"
        self.plink_fn = os.path.join(self.args.output_dir, plink_fn)
        self.plink_commands_f = None
            
    def writeAnalysisReport(self):
        report_fn = os.path.join(self.args.output_dir, self.args.report_fn)
        with open(report_fn, 'w') as out_f:
            out_f.write("Program version: %s\n" % __version__)
            out_f.write("Running date: %s\n\n" % time.ctime())
            out_f.write("Input dataset:\t%s\n" % self.args.input)
            #no_removed = self.failed_samples.getFailCount("sex")
            #out_f.write("Removed %d samples in sex check\n" % (no_removed))
            out_f.write("Missingness rate:\t%3.5f\n" % self.args.missingness)
            out_f.write("Heterozygosity rate:\t%3.5f * SD\n" % self.args.heterozygosity)
            #no_removed = self.failed_samples.getFailCount("miss_and_het")
            #out_f.write("Removed %d samples in missingness and heterozygosity check\n" % (no_removed))
            out_f.write("Pi-hat value:\t%3.5f\n" % self.args.pi_hat)
            #no_removed = self.failed_samples.getFailCount("ibd")
            #out_f.write("Removed %d samples in contamination check\n" % (no_removed))
            out_f.write("Principal component 1 cutoff:\t%3.5f * SD\n" % self.args.C1)
            out_f.write("Principal component 2 cutoff:\t%3.5f * SD\n" % self.args.C2)
            out_f.write("Maximum number of MDS rounds: %d\n" % self.args.max_mds_rounds)
            out_f.write("Principal component 1 cutoff string: %s\n" % self.args.C1_string)
            out_f.write("Principal component 2 cutoff string: %s\n" % self.args.C2_string)
            #no_removed = self.failed_samples.getFailCount("mds")
            #out_f.write("Removed %d samples in multidimensional scaling check\n" % (no_removed))
            out_f.write("Call rate limit:\t%3.5f\n" % self.args.geno)
            out_f.write("HWE limit:\t%3.5f\n" % self.args.hwe)
            out_f.write("MAC limit for imputation:\t%3.5f\n" % self.args.mac)
    
    
    def main(self):
        self.plink_commands_f = open(self.plink_fn, 'w')
        # filter sex
        if self.args.noX:
            sex_filtered = self.args.input
        else:
            filter = FilterBySex(self.args, plink_commands_f = self.plink_commands_f)
            sex_filtered = filter.runComponent(self.args.input)
        
        # filter missingness / heterozygosity
        filter = FilterHeterozygosityAndMissingRate(self.args, plink_commands_f = self.plink_commands_f)
        het_missing_filtered = filter.runComponent(sex_filtered)
        
        # IBS calculation and contamination filtering
        ibs_filter = IBSCalculation(self.args, plink_commands_f = self.plink_commands_f)
        before_mds_fn, pihat_fn = ibs_filter.runComponent(het_missing_filtered, filterContaminations = True)
        
        # remove related for MDS calculations
        remove_related = RemoveRelatedFromMDS(self.args, plink_commands_f = self.plink_commands_f)
        mds_fn = remove_related.runComponent(before_mds_fn, pihat_fn)
        
        # loop the MDS
        mds_filter = FilterMDS(self.args, plink_commands_f = self.plink_commands_f)
        round_no = 0
        final = False
        while not final:
            round_no += 1
            genome_fn = ibs_filter.runComponent(mds_fn, round_no = round_no)
            mds_fn, no_overlimit = mds_filter.runComponent(mds_fn, genome_fn, round_no)
            # redo genome for next round of MDS
            if (no_overlimit == 0) or (round_no > self.args.max_mds_rounds):
                final = True

        # QC snp exclusion
        filter_snps = FilterSNPExclusion(self.args, plink_commands_f = self.plink_commands_f)
        qc_done_fn = filter_snps.runComponent(before_mds_fn)
        
        # prepare dataset for imputation
        impu_filter = FilterForImputation(self.args, plink_commands_f = self.plink_commands_f)
        for_imputation_fn = impu_filter.runComponent(qc_done_fn)
        # result report
        self.writeAnalysisReport()
        self.plink_commands_f.close()
       
        
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
    #parser.add_argument("--pi_hat_removal_limit", help = "Pi-hat limit for removing samples", type = float, default = 0.95)
    parser.add_argument("--hwe", help = "HWE limit", type = float, default = 0.000001)
    parser.add_argument("--geno", help = "Call rate limit", type = float, default = 0.02)
    parser.add_argument("--mac", help = "MAC treshold for imputation", type = float, default = 3)
    parser.add_argument("--repeat", help = "Removing duplicates who's entry is repeated more than this limit (due to contamination).",
                        type = int, default = 15)
    parser.add_argument("--noX", help = "No X chromosome in the data - skip the sex checks", default = False, 
                        action = "store_true")
#    parser.add_argument("--plink_path", help = "Path to the plink executable", default = "/apps/genetics/bin/plink-1.9")
    parser.add_argument("--report_fn", help = "Analysis report file name", default = "analysis_report.txt")
    parser.add_argument("--high_ld_regions", help = "Path to the file used in duplicate filtering",
                        default = "/fs/projects/finngen/misc/high-LD-regions.txt")
    parser.add_argument("--version", help = "Shows the version number of the program",
                        action = 'version', version = "%(prog)s {version}".format(version=__version__))
    return parser.parse_args(args)

if __name__=='__main__':
    app = GWAS_QC()
    app.main()
    