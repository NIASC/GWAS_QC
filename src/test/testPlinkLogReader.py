'''
Created on 8.11.2016

@author: paakkone
'''
import unittest
from core.plinkLogReader import PlinkLogReader
import os

class Test(unittest.TestCase):
    def setUp(self):
        self.log_fns = ("check_sex.log",
                        "removed_sex.log",
                        "check_missingness.log",
                        "removed_missingness.log",
                        "check_heterozygosity.log",
                        "removed_heterozygosity.log",
                        "removed_heterozygosity_check_pruning.log",
                        "removed_heterozygosity_check_genome.log",
                        "removed_ibs.log",
                        "removed_mds_rel.log",
                        "removed_mds_rel_check_pruning1.log",
                        "removed_mds_rel_check_genome1.log",
                        "calc_mds_round_1.log",
                        "removed_mds_final.log",
                        "removed_all_mds.log",
                        "SUMMIT_CR_98_gender_monomorphic_missing_HWE_HET_REL_FINAL_TOP_b37_QCed.log",
                        "calc_freq.log",
                        "SUMMIT_CR_98_gender_monomorphic_missing_HWE_HET_REL_FINAL_TOP_b37_for_imputation.log",
                        "summit_denials_removed_for_imputation.log")


    def tearDown(self):
        pass


    def testName(self):
        for fn in self.log_fns:
            fn = os.path.join("test_data", fn)
            print ("FN:",fn)
            log = PlinkLogReader(fn)
            self.assertNotEqual(log.plink_ver, None)
            self.assertNotEqual(log.start_variants, None)
            self.assertNotEqual(log.start_samples, None)
            #print ("start_variants: %d\tstart_samples: %d" % (log.start_variants, log.start_samples) )
            if fn.find("missing") == -1 and fn.find("calc_freq") == -1:
                self.assertNotEqual(log.end_variants, None)
                self.assertNotEqual(log.end_samples, None)
                #print ("end_variants: %d\t end_samples: %d" % (log.end_variants, log.end_samples))
                self.assertTrue(log.start_variants >= log.end_variants)
                self.assertTrue(log.start_samples >= log.end_samples)
            print ("Removed %d samples and %d variants" % (log.getNoRemovedSamples(), log.getNoRemovedVariants()))
            if log.hwe_removed != None:
                print ("fn: %s\t HWE removed: %d" % (fn, log.hwe_removed))
        
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()