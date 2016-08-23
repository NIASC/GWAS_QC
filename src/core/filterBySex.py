'''
Created on 19.8.2016

@author: paakkone
'''

import os
from component import Component

class FilterBySex(Component):

    def __init__(self, args, failed_fn = None, plink_commands_f = None):
        self.fail_type = "sex"
        super().__init__(args, failed_fn = failed_fn, plink_commands_f = plink_commands_f)
        
    def runPlink(self, dataset_fn):
        #root, head = os.path.split(dataset_fn)
        out_fn = "check_%s" % self.fail_type
        out_fn = os.path.join(self.args.output_dir, out_fn)
        switches = ["--check-sex"]
        self.plinkRunner.runPlinkCommand(switches, dataset_fn, out_fn)
        return out_fn
    
    def findFailedSamples(self, out_fn):
        fn = "%s.sexcheck" % (out_fn)
        with open(fn) as in_f:
            for line in in_f:
                words = line.split()
                if words[4] == "PROBLEM":
                    sex1 = words[2]
                    sex2 = words[3]
                    if sex1 in ["1", "2"] and sex2 in ["1", "2"]:
                        fam_id = words[0]
                        p_id = words[1]
                        self.failed_samples.append([fam_id, p_id])
                        #self.failed_samples.addSample(fam_id, p_id, self.fail_type)
    
    def runComponent(self, dataset_fn):                    
        # prepare switches and file names
        print ("Run Component: %s" % dataset_fn)
        out_fn = self.runPlink(dataset_fn)
        # filter for failures
        self.findFailedSamples(out_fn)
        # write out failures
        self.writeFailedSamples()
        # remove samples and write new file set
        new_dataset_fn = self.removeIndividuals(dataset_fn)
        # return new file set name
        return new_dataset_fn
    
    