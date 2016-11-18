'''
Created on 18.8.2016

@author: paakkone
'''

import os
from plotter import Plotter
from collections import OrderedDict
from component import Component
from utils import getUpperAndLowerLimits

class FilterHeterozygosityAndMissingRate(Component):
    
    def __init__(self, args, failed_fn = None, plink_commands_f = None):
        self.fail_type = "miss_and_het"
        super().__init__(args, failed_fn = failed_fn, plink_commands_f = plink_commands_f)
        self.imiss_values = OrderedDict()
        self.het_values = OrderedDict()
        # geno and hwe switch values used for heterozygosity calculation
        self.geno = 0.02
        self.hwe = 0.000001
        
    def readIMiss(self, fn):
        """
        Reads IMiss file, writes entries which exceed the missingness limit to
        remove_imiss file and store rest to self.imiss_values ordered dictionary.
        Keys of the dictionary are combined fid & iid values.
        """
        print ("readIMiss fn: %s" % fn)
        root, head = os.path.split(fn)
        imiss_fn = fn + ".imiss"
        remove_fn = "remove_imiss.txt"
        remove_fn = os.path.join(self.args.output_dir, remove_fn)
        remove_f = open(remove_fn, 'w')       
        with open(imiss_fn) as in_f:
            for line in in_f:
                words = line.split()
                if words[0] != "FID":
                    imiss_val = words[5]
                    if float(imiss_val) >= self.args.missingness:
                        # failed missingness
                        remove_f.write("%s\t%s\n" % (words[0], words[1]))
                    else:               
                        key = words[0] + "&" + words[1]
                        self.imiss_values[key] = imiss_val
        remove_f.close()
        return remove_fn
                    
                    
    def runPlink(self, dataset_fn):
        """
        Missingness calculations
        """
        root, head = os.path.split(dataset_fn)
        out_fn = "check_missingness"
        out_fn = os.path.join(self.args.output_dir, out_fn)
        switches = ["--missing","--allow-no-sex"]
        # calculate missingness
        self.plinkRunner.runPlinkCommand(switches, dataset_fn, out_fn)
        self.tempFilesCreated(out_fn)
        return out_fn
    
    def runPlink2(self, imiss_removed_fn):
        """
        Heterozygosity calculations
        """
        root, head = os.path.split(imiss_removed_fn)
        out_fn = "check_heterozygosity"
        out_fn = os.path.join(self.args.output_dir, out_fn)
        switch1 = "--geno %s" % self.geno
        switch2 = "--hwe %s" % self.hwe
        switches = [switch1, switch2, "--het","--allow-no-sex"]
        # calculate missingness
        self.plinkRunner.runPlinkCommand(switches, imiss_removed_fn, out_fn)
        self.tempFilesCreated(out_fn)
        return out_fn
        
    def findFailedSamples(self, dataset_fn, remove_het_f):
        fn = "%s.het" % (dataset_fn)
        # delete values which have non_miss_no = 0 from both heterozygosity and missingness values
        with open(fn) as in_f:
            for line in in_f:
                words = line.split()
                if words[0] != "FID":
                    fam_id = words[0]
                    p_id = words[1]
                    key = fam_id + "&" + p_id
                    hom = float(words[2])
                    non_miss_no = float(words[4])
                    if non_miss_no == 0:
                        #indexes_to_delete.append(i-1)
                        self.failed_samples.append([fam_id, p_id])
                        #self.failed_samples.addSample(fam_id, p_id)
                        # delete this one also from imiss values
                        del self.imiss_values[key]
                        # write the deleted sample id to removed file
                        remove_het_f.write("%s\t%s\n" % (fam_id, p_id))
                    else:
                        obs_het = (non_miss_no - hom) / (non_miss_no * 1.0)
                        self.het_values[key] = obs_het

    def findValuesGreaterThanHeterozygosityLimit(self, remove_het_f):
        std_limit = float(self.args.heterozygosity)
        het_vals = [self.het_values[key] for key in self.het_values.keys()]
        # color map
        cm = ["green"] * len(het_vals)
        lower, upper = getUpperAndLowerLimits(het_vals, std_limit)
        i = 0
        for key in self.het_values.keys():
            het_val = self.het_values[key]
            if het_val > upper or het_val < lower:
                fam_id, p_id = key.split("&")
                remove_het_f.write("%s\t%s\n" % (fam_id, p_id))
                #self.failed_samples.addSample(fam_id, p_id, fail_type)
                cm[i] = "red"
            i += 1
        return cm
    
    def plot(self, color_map):
        out_fn = os.path.join(self.args.output_dir, "missingness_and_heterozygosity.png")
        std = float(self.args.heterozygosity)
        miss_limit = float(self.args.missingness)
        het_vals = [float(self.het_values[key]) for key in self.het_values.keys()]
        miss_vals = [float(self.imiss_values[key]) for key in self.imiss_values.keys()]
        figure_title = "Heterozygosity rate and missingness"
        x_title = "Proportion of missing genotypes"
        y_title = "Heterozygosity rate"
        plotter = Plotter()
        plotter.plotScatterPlot(out_fn, miss_vals, het_vals, x_val_limit = miss_limit, 
                                  y_std_limit = std, figure_title = figure_title, 
                                  x_title = x_title, y_title = y_title, x_log_scale = True, col_map = color_map)

    def runComponent(self, dataset_fn):
        # do plink calculations
        out_fn = self.runPlink(dataset_fn)
        # read the imiss file
        remove_imiss_fn = self.readIMiss(out_fn)
        # remove bad missingness values and create new dataset
        #print ("Removing imiss failures")
        imiss_removed_fn = self.removeIndividuals(dataset_fn, failed_samples = remove_imiss_fn, output_fn = "removed_missingness")
        #print ("Removed, %s" % imiss_removed_fn)
        # do second calculation
        het_fn = self.runPlink2(imiss_removed_fn)
        #self.calculate(switches, imiss_removed_fn, out_fn)
        remove_het_fn = os.path.join(self.args.output_dir, "remove_het.txt")
        remove_het_f = open(remove_het_fn, 'w')
        # filter for failures
        self.findFailedSamples(het_fn, remove_het_f)
        # find samples over heterozygosity limit
        cm = self.findValuesGreaterThanHeterozygosityLimit(remove_het_f)
        # remove samples and write new file set
        remove_het_f.close()
        new_dataset_fn = self.removeIndividuals(imiss_removed_fn, failed_samples = remove_het_fn, output_fn = "removed_heterozygosity")
        # plot missingness_heterozygosity plot
        self.plot(cm)
        # read log file for reporting of number of samples
        # there are two log files in this..
        self.log.readLogFile(new_dataset_fn + ".log")
        # return new file set name
        return new_dataset_fn
