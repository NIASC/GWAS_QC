'''
Created on 19.8.2016

@author: paakkone
'''

import os
from component import Component
from plotter import Plotter
from utils import getUpperAndLowerLimits

class FilterMDS(Component):
    
    def __init__(self, args, failed_fn = None, plink_commands_f = None):
        self.fail_type = "mds"
        super().__init__(args, failed_fn = failed_fn, plink_commands_f = plink_commands_f)
  
    def runPlink(self, dataset_fn, gen_fn, round_no):
        gen_fn = gen_fn + ".genome"
        switch1 = "--read-genome %s" % gen_fn
#         if round_no == 1:
#             output_fn = dataset_fn + "_%s" % self.fail_type
#         else:
#             # remove _idb
#             last_dash = dataset_fn.rfind("_")
#             # remove previous rounds _mdsrun#
#             last_dash = dataset_fn[:last_dash].rfind("_")
#             output_fn = dataset_fn[:last_dash] + "_%s" % self.fail_type
        output_fn = "calc_mds_round_%d" % round_no
        output_fn = os.path.join(self.args.output_dir, output_fn)
        self.plinkRunner.runPlinkCommand([switch1, "--cluster", "--mds-plot 10 "], dataset_fn, output_fn)
        return output_fn
        
    def runComponent(self, dataset_fn, gen_fn, round_no):
        # run plink
        ds_fn = self.runPlink(dataset_fn, gen_fn, round_no)
        # filter for failures
        c1_vals, c2_vals, c1_limit, c2_limit, cm, failed_no = self.findFailedSamples(ds_fn, round_no)
        # write out failures
        if failed_no > 0:
            fail_fn = "remove_mds_round%d.txt" % round_no
            out_fn = "removed_mds_round%d" % round_no
        else:
            fail_fn = "remove_mds_all.txt"
            out_fn = "removed_mds_final"
        fail_fn = os.path.join(self.args.output_dir, fail_fn)
        out_fn = os.path.join(self.args.output_dir, out_fn)
        self.writeFailedSamples(fail_fn)
        # remove samples and write new file set
        new_dataset_fn = self.removeIndividuals(dataset_fn, failed_samples = fail_fn, output_fn = out_fn)
        # plotting
        self.plot(c1_vals, c2_vals, c1_limit, c2_limit, round_no, cm)
        # return new file set name
        return new_dataset_fn, failed_no
    
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
    
    
    def findFailedSamples(self, dataset_fn, round_no):
        """
        Reads MDS file and finds values violating C1 and C2 limits. Returns values and color map for plotting.
        """
        mds_fn = dataset_fn + ".mds"
        # read MDS file
        failed_no = 0
        header = True
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
                    mds_data.append([fid, iid, c1, c2])
        
        # calculate limits
        c1_vals = [c1 for fid, iid, c1, c2 in mds_data]
        c2_vals = [c2 for fid, iid, c1, c2 in mds_data]
        c1_limit, c2_limit = self.getMSDLimit(round_no)
        c1_std_lower_limit, c1_std_upper_limit = getUpperAndLowerLimits(c1_vals, c1_limit)
        c2_std_lower_limit, c2_std_upper_limit = getUpperAndLowerLimits(c2_vals, c2_limit)
        
        # find violations   
        col_map = ["green"] * len(mds_data)     
        i = 0
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
                self.failed_samples.append([fid, iid])
#                self.failed_samples.addSample(fid, iid)
                col_map[i] = "red"
                failed_no += 1
            i += 1
        return c1_vals, c2_vals, c1_limit, c2_limit, col_map, failed_no
    
    def plot(self, c1_vals, c2_vals, c1_limit, c2_limit, round_no, cm):
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
        plotter.plotScatterPlot(out_fn, c1_vals, c2_vals, x_std_limit = c1_limit, y_std_limit = c2_limit, 
                                  figure_title = figure_title, x_title = x_title,
                                  y_title = y_title, plot_std_y_axis = True, col_map = cm)
    
    