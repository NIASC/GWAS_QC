'''
Created on 15.12.2015

@author: paakkone
'''
import unittest
from core.gwas_qc import GWAS_QC

class Test(unittest.TestCase):

    def setUp(self):
        self.dataset_name = "test_data/raw-GWA-data"
        self.output_dir = "test_data"
        self.input_params = ["--input", self.dataset_name, "--output_dir", self.output_dir, "--plink_path", ""]
        self.gwas_qc = GWAS_QC(self.input_params)

    def load_dataset2(self, ds_name = None):
        if ds_name == None:
            ds_name = "C:/Users/paakkone/test_data/qc_filtering/original/PLUSSTRAND_GeneRISK_DHR_MRPred"
        output_dir = "C:/Users/paakkone/test_data/qc_filtering/output/"
        input_par = ["--input", ds_name, "--output_dir", output_dir, "--plink_path", "C:/apps/plink-1_90/plink.exe"]
        return GWAS_QC(input_par)

    

#     def testScatterPlot(self):
#         self.gwas_qc.filterHeterozygosityAndMissingRate()

#     def testFilterBySex(self):
#         qc = self.load_dataset2()
#         qc.filterBySex()
        
#     def testFilterByMissingnessAndHeterotsygosity(self):
#         qc = self.load_dataset2("C:/Users/paakkone/test_data/qc_filtering/output/missing2")
#         qc.filterHeterozygosityAndMissingRate()
        
    def testMain(self):
        qc = self.load_dataset2()
        qc.main()
#    def testHistogramPlot(self):
#        self.gwas_qc.filterMissingDataRate()
        
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()