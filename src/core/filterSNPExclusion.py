'''
Created on 22.8.2016

@author: paakkone
'''

import os
from component import Component

class FilterSNPExclusion(Component):
    '''
    classdocs
    '''

    def __init__(self, args, failed_fn = None, plink_commands_f = None):
        self.fail_type = "snp"
        super().__init__(args, failed_fn = failed_fn, plink_commands_f = plink_commands_f)
        
    def runPlink(self, input_fn):
        """
        SNP exclusion
        """
        orig_input_fn = self.args.input
        in_path, in_fn = os.path.split(self.args.input)
        out_fn = in_fn + "_QCed"
        out_fn = os.path.join(self.args.output_dir, out_fn)
        switch1 = "--geno %s" % self.args.geno
        switch2 = "--hwe %s" % self.args.hwe
        self.plinkRunner.runPlinkCommand([switch1, switch2, "--make-bed"], input_fn, out_fn)
        # don't store these files to files created to avoid them being deleted as temporary ones
        return out_fn
   
   
    def runComponent(self, ds_name):
        # remove samples found in mds runs
        remove_fn = "remove_mds_all.txt"
        remove_fn = os.path.join(self.args.output_dir, remove_fn)
        remove_ds_fn = "removed_all_mds"
        remove_ds_fn = self.removeIndividuals(ds_name, failed_samples = remove_fn, output_fn = remove_ds_fn)
        # qc snp filtering
        qc_done_fn = self.runPlink(remove_ds_fn)
        # read log file for reporting of number of samples
        self.log.readLogFile(qc_done_fn + ".log")
        return qc_done_fn
        
    def findFailedSamples(self):
        pass
    
#