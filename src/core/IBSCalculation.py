'''
Created on 18.8.2016

@author: paakkone
'''

import os
from component import Component

class IBSCalculation(Component):
    def __init__(self, args, failed_fn = None, plink_commands_f = None):
        self.fail_type = "ibs"
        super().__init__(args, failed_fn = failed_fn, plink_commands_f = plink_commands_f)
        self.geno = 0.02
        self.hwe = 0.000001
        
    def runPlink(self, in_fn, round_no):
        """
        Pruning calculations
        """
        switch1 = "--geno %s" % self.geno
        switch2 = "--hwe %s" % self.hwe
        switch3 = "--exclude %s --range" % self.args.high_ld_regions
        out_fn = in_fn + "_check_pruning"
        #out_fn = "check_pruning"
        if round_no != None:
            out_fn += "%d" % round_no
        self.plinkRunner.runPlinkCommand([switch1,switch2, switch3,"--indep-pairwise 50 5 0.2"], in_fn, out_fn)
        return out_fn
        
    def runPlink2(self, in_fn, pruned_fn, round_no):
        """
        genome calculations
        """
        extract_str = "--extract %s.prune.in" % pruned_fn
        out_fn = in_fn + "_check_genome"
        #out_fn = "check_genome"
        if round_no != None:
            out_fn += "%d" % round_no
        self.plinkRunner.runPlinkCommand([extract_str, "--genome"], in_fn, out_fn)
        return out_fn
    
    def runComponent(self, in_fn, filterContaminations = False, round_no = None):
        # pruning calculations
        pruned_fn = self.runPlink(in_fn, round_no)
        # genotype calculations
        geno_fn = self.runPlink2(in_fn, pruned_fn, round_no)
        if filterContaminations == False:
            print ("Returning IBS calc run: geno_fn %s" % (geno_fn))
            return geno_fn
        pihat_fn = self.findFailedSamples(geno_fn)
        # write out failures
        self.writeFailedSamples()
        # remove samples and write new file set
        new_dataset_fn = self.removeIndividuals(in_fn)
        # return new file set name
        print ("Returning IBS calc run: new_dataset_fn %s  pihat_fn %s" % (new_dataset_fn, pihat_fn))
        return new_dataset_fn, pihat_fn
    
    def findFailedSamples(self, dataset_fn):
        """
        Finds contaminated samples (samples which repeat more than given limit) and removes them.
        Lists IDs which have pi_hat values larger than given limit into a file.
        """
        delim = "&"
        remove = []
        key_counter = {}
        repeat_limit = int(self.args.repeat)
        pihat_over_limit = []
        header_row = True
        genome_fn = "%s.genome" % dataset_fn
        pihat_fn = "PI_hat_over_%1.2f.txt" % self.args.pi_hat
        pihat_fn = os.path.join(self.args.output_dir, pihat_fn)
        pihat_f = open(pihat_fn, 'w')
        with open(genome_fn) as in_f:
            for line in in_f:
                words = line.split()
                if header_row:
                    # write header to pihat_fn
                    header_row = False
                    pihat_f.write("%s\t%s\t%s\t%s\t%s\n" % (words[0],words[1],words[2],words[3],words[9].strip()))
                else:
                    fid1 = words[0]
                    iid1 = words[1]
                    fid2 = words[2]
                    iid2 = words[3]
                    key1 = fid1 + delim + iid1
                    key2 = fid2 + delim + iid2
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
                                    fid, pid = key.split(delim)
                                    self.failed_samples.append([fid, pid])
                                    #self.failed_samples.addSample(fid, pid)
        # go through the pihat list again and write out those over pihat limit
        # cannot do already on previous round as the keys to be removed is not known
        # before whole file has been checked
        for key1, key2, pi_hat in pihat_over_limit:
            if (key1 in remove) or (key2 in remove):
                pass
            else:
                fid1, pid1 = key1.split(delim)
                fid2, pid2 = key2.split(delim)
                pihat_f.write("%s\t%s\t%s\t%s\t%s\n" % (fid1, pid1, fid2, pid2, pi_hat))
        pihat_f.close()
        return pihat_fn

