'''
Created on 22.8.2016

@author: paakkone
'''
import os
from collections import OrderedDict
from component import Component

class FilterForImputation(Component):
    '''
    classdocs
    '''

    def __init__(self, args, failed_fn = None, plink_commands_f = None):
        self.fail_type = "impu"
        super().__init__(args, failed_fn = failed_fn, plink_commands_f = plink_commands_f)
        
    def runPlink(self, input_fn):
        """
        Calculate frequencies
        """
        freq_fn = "calc_freq"
        freq_fn = os.path.join(self.args.output_dir, freq_fn)
        self.plinkRunner.runPlinkCommand(["--freq"], input_fn, freq_fn)
        return freq_fn
     
    def runPlink2(self, input_fn, fail_fn):
          
        orig_input_fn = self.args.input
        in_path, in_fn = os.path.split(self.args.input)
        out_fn = in_fn + "_for_imputation"
        out_fn = os.path.join(self.args.output_dir, out_fn)
          
        switch1 = "--exclude %s" % fail_fn
        self.plinkRunner.runPlinkCommand([switch1, "--make-bed"], input_fn, out_fn)
        return out_fn
    
    
    def runComponent(self, ds_name):
        freq_fn = self.runPlink(ds_name)
        
        #return qc_done_fn
        failed_fn = self.findFailedSamples(ds_name, freq_fn)
        
        # exclude failing snips
        impu_fn = self.runPlink2(ds_name, failed_fn)
        return impu_fn
        
    def findFailedSamples(self, input_fn, freq_fn):
        """
        Use bim and freq file to filter non-autosomal snips and snps with 
        MAC value less than  treshold
        """
        #freq_fn = freq_fn + "_freq"
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
        return fail_fn
    
    
    