'''
Created on 8.11.2016

@author: paakkone
'''

class PlinkLogReader(object):
    '''
    Class to read Plink log files and extract data from it
    '''
    def __init__(self):
        '''
        Constructor
        '''
        self.plink_ver = None
        self.start_variants = None
        self.start_samples = None
        self.end_variants = None
        self.end_samples = None
        self.hwe_removed = None
        
    def readLogFile(self, log_fn):
        #print ("Reading log %s" % log_fn)
        with open(log_fn) as log_f:
            for line in log_f:
                line = line.strip()
                words = line.split()
                if line.startswith("PLINK"):
                    self.plink_ver = line.strip()
                if line.endswith("variants loaded from .bim file."):
                    self.start_variants = int(words[0])
                if line.endswith("loaded from .fam.") and line.find("people") != -1:
                    self.start_samples = int(words[0])
                if line.endswith("pass filters and QC."):
                    self.end_variants = int(words[0])
                    self.end_samples = int(words[3])
                if line.endswith("variants removed due to Hardy-Weinberg exact test."):
                    self.hwe_removed = int(words[1])
        #print ("start samples: %s \nstart variants %s" % (self.start_samples, self.start_variants))
        #print ("end samples: %s \nend variants %s" % (self.end_samples, self.end_variants))
        
    def getNoStartSamples(self):
        return self.start_samples
    
    def getNoStartVariants(self):
        return self.start_variants
        
    def getNoRemovedSamples(self):
        if self.end_samples == None:
            return 0
        return self.start_samples - self.end_samples
    
    def getNoRemovedVariants(self):
        if self.end_variants == None:
            return 0
        return self.start_variants - self.end_variants
        
    def getNoHWERemoved(self):
        if self.hwe_removed == None:
            return 0
        return self.hwe_removed
                
    def getNoEndVariants(self):
        if self.end_variants == None:
            return self.start_variants
        return self.end_variants
    
    def getNoEndSamples(self):
        if self.end_samples == None:
            return self.start_samples
        return self.end_samples
    