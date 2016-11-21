#!/apps/python3.5/anaconda3/bin/python
'''
Created on 18.8.2016

@author: paakkone
'''
__version__= "1.3"

import sys
import os
import time
import datetime
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
        self.start_time = datetime.datetime.now()
        self.args = parseArguments(args)
        if not os.path.exists(self.args.output_dir):
            os.makedirs(self.args.output_dir)
        plink_fn = "plink_commands.txt"
        self.plink_fn = os.path.join(self.args.output_dir, plink_fn)
        self.plink_commands_f = None
        self.samples = {}
        self.variants = {}
        self.files_to_delete = []
        
    def writeAnalysisReport(self):
        report_fn = os.path.join(self.args.output_dir, self.args.report_fn)
        with open(report_fn, 'w') as out_f:
            out_f.write("Program version: %s\n" % __version__)
            out_f.write("Running date: %s\n\n" % time.ctime())
            out_f.write("Working directory: \n%s\n" % os.getcwd())
            out_f.write("Commandline:\n%s\n\n" % " ".join(sys.argv))
            self.printParameters(out_f)
            self.printRemoved(out_f)
            out_f.write("%s\n" % (self.getElapsedTime()))
    
    def printParameters(self, out_f):
        out_f.write("PARAMETERS\n")
        out_f.write("Input dataset:\t%s\n" % self.args.input)
        out_f.write("Missingness rate:\t%3.5f\n" % self.args.missingness)
        out_f.write("Heterozygosity rate:\t%3.5f * SD\n" % self.args.heterozygosity)
        out_f.write("Pi-hat value:\t%3.5f\n" % self.args.pi_hat)
        out_f.write("Principal component 1 cutoff:\t%3.5f * SD\n" % self.args.C1)
        out_f.write("Principal component 2 cutoff:\t%3.5f * SD\n" % self.args.C2)
        out_f.write("Maximum number of MDS rounds: %d\n" % self.args.max_mds_rounds)
        out_f.write("Principal component 1 cutoff string: %s\n" % self.args.C1_string)
        out_f.write("Principal component 2 cutoff string: %s\n" % self.args.C2_string)
        out_f.write("Call rate limit:\t%3.5f\n" % self.args.geno)
        out_f.write("HWE limit:\t%3.8f\n" % self.args.hwe)
        out_f.write("MAC limit for imputation:\t%3.5f\n" % self.args.mac)
    
    def printRemoved(self, out_f):
        out_f.write("\nCALCULATIONS\n")
        out_f.write("Number of samples in the dataset: %s\n" % self.samples["start"])
        out_f.write("Number of variants in the dataset: %s\n" % self.variants["start"])
        if "sex" in self.samples.keys():
            out_f.write("Removed %s samples in sex check\n" % (self.samples["sex"]))
        out_f.write("Removed %s samples in heterozygosity and missingness checks\n" % self.samples["hetmiss"])
        out_f.write("Removed %s samples in IBS calculations\n" % self.samples["ibs"])
        out_f.write("Removed %s samples in MDS calculations\n" % self.samples["mds"])
        #out_f.write("Removed %s variants in MDS calculations\n" % self.variants["mds"])
        qc_samples_removed = self.samples["start"] - self.samples["qc"]
        qc_variants_removed = self.variants["start"] - self.variants["qc"]
        out_f.write("Number of samples in the QCd dataset: %s (removed %s)\n" % (self.samples["qc"], qc_samples_removed))
        out_f.write("Number of variants in the QCd dataset: %s (removed %s)\n" % (self.variants["qc"], qc_variants_removed))
        impu_samples_removed = self.samples["start"] - self.samples["impu"]
        impu_variants_removed = self.variants["start"] - self.variants["impu"]
        out_f.write("Number of samples in the dataset for imputation: %s (removed %s)\n" % (self.samples["impu"], impu_samples_removed))
        out_f.write("Number of variants in the dataset for imputation: %s (removed %s)\n" % (self.variants["impu"], impu_variants_removed))
    
        
    def getElapsedTime(self):
        end_time = datetime.datetime.now()
        elapsed_time = end_time - self.start_time
        hours = int(elapsed_time / datetime.timedelta(hours=1))
        minutes = int(elapsed_time / datetime.timedelta(minutes = 1)) % 60
        seconds = int(elapsed_time.total_seconds() %  60)
        time_str = "\nTime elapsed: %d hours, %d minutes, %d seconds" % (hours, minutes, seconds)
        return time_str
    
    def deleteIntermediateFiles(self):
        """
        Deletes the files created during calculation but not necessary in the end
        Keeps the QC and for_imputation datasets
        """
        uniq_files = set(self.files_to_delete)
        print ("Deleting %d intermediate files" % len(uniq_files))
        for fn in uniq_files:
            #print ("Deleting file %s" % fn)
            os.remove(fn)
    
    def main(self):
        self.plink_commands_f = open(self.plink_fn, 'w')
        # filter sex
        if self.args.noX:
            sex_filtered = self.args.input
        else:
            sex_filter = FilterBySex(self.args, plink_commands_f = self.plink_commands_f)
            sex_filtered = sex_filter.runComponent(self.args.input)
            self.files_to_delete.extend(sex_filter.getTempFiles())
            self.samples["start"] = sex_filter.getNoStartSamples()
            self.variants["start"] = sex_filter.getNoStartVariants()
            self.samples["sex"] = sex_filter.getNoRemovedSamples()
        # filter missingness / heterozygosity
        hetmiss_filter = FilterHeterozygosityAndMissingRate(self.args, plink_commands_f = self.plink_commands_f)
        het_missing_filtered = hetmiss_filter.runComponent(sex_filtered)
        self.files_to_delete.extend(hetmiss_filter.getTempFiles())
        if self.args.noX:
            self.samples["start"] = hetmiss_filter.getNoStartSamples()
            self.variants["start"] = hetmiss_filter.getNoStartVariants()
        self.samples["hetmiss"] = sex_filter.getNoEndSamples() - hetmiss_filter.getNoEndSamples()
        # IBS calculation and contamination filtering
        ibs_filter = IBSCalculation(self.args, plink_commands_f = self.plink_commands_f)
        before_mds_fn, pihat_fn = ibs_filter.runComponent(het_missing_filtered, filterContaminations = True)
        self.files_to_delete.extend(ibs_filter.getTempFiles())
        self.samples["ibs"] = ibs_filter.getNoRemovedSamples()
        
        # remove related for MDS calculations
        remove_related = RemoveRelatedFromMDS(self.args, plink_commands_f = self.plink_commands_f)
        mds_fn = remove_related.runComponent(before_mds_fn, pihat_fn)
        self.files_to_delete.extend(remove_related.getTempFiles())
        
        # loop the MDS
        mds_filter = FilterMDS(self.args, plink_commands_f = self.plink_commands_f)
        round_no = 0
        final = False
        while not final:
            round_no += 1
            genome_fn = ibs_filter.runComponent(mds_fn, round_no = round_no)
            self.files_to_delete.extend(ibs_filter.getTempFiles())
            mds_fn, no_overlimit = mds_filter.runComponent(mds_fn, genome_fn, round_no)
            self.files_to_delete.extend(mds_filter.getTempFiles())
            # redo genome for next round of MDS
            if (no_overlimit == 0) or (round_no > self.args.max_mds_rounds):
                final = True
        self.samples["mds"] = mds_filter.getNoFinalMDSRemovals()
        #self.variants["mds"] = ibs_filter.getNoEndVariants() - mds_filter.getNoEndVariants()
        # QC snp exclusion
        filter_snps = FilterSNPExclusion(self.args, plink_commands_f = self.plink_commands_f)
        qc_done_fn = filter_snps.runComponent(before_mds_fn)
        self.files_to_delete.extend(filter_snps.getTempFiles())
        self.samples["qc"] = filter_snps.getNoEndSamples()
        self.variants["qc"] = filter_snps.getNoEndVariants()
        
        # prepare dataset for imputation
        impu_filter = FilterForImputation(self.args, plink_commands_f = self.plink_commands_f)
        for_imputation_fn = impu_filter.runComponent(qc_done_fn)
        self.samples["impu"] = impu_filter.getNoEndSamples()
        self.variants["impu"] = impu_filter.getNoEndVariants()
        # result report
        self.writeAnalysisReport()
        self.plink_commands_f.close()
        if self.args.del_intermediate:
            self.deleteIntermediateFiles()
        
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
    parser.add_argument("--del_intermediate", help = "Delete intermediate files", default = False,
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
    