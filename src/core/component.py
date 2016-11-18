'''
Created on 15.8.2016

@author: paakkone
'''

import os
import abc
import glob
from plotter import Plotter
from utils import PlinkRunner
from plinkLogReader import PlinkLogReader

class Component(object, metaclass=abc.ABCMeta):
    '''
    General component object from which actual components inherit from
    '''
    def __init__(self, args, failed_fn = None, plink_commands_f = None):
        '''
        Takes the arguments (self.args) of the main program
        '''
        self.args = args
        #self.fail_type = fail_type
        self.failed_samples = []
        self.failed_fn = failed_fn
        if failed_fn == None:
            remove_fn = "remove_%s.txt" % self.fail_type
            self.failed_fn = os.path.join(self.args.output_dir, remove_fn)
            #self.failed_fn = "removed_%s.txt" % self.fail_type
        self.plinkRunner = PlinkRunner(plink_commands_f)
        # log file will be read at the end of run component
        self.log = PlinkLogReader()
        # keep track of files created so that they can be deleted if intermediate files take too much space
        # don't include final results, removal lists or plots
        self.temp_files = []
   
    @abc.abstractmethod
    def runPlink(self):
        raise NotImplementedError("Each class needs to implement runPlink")
   
    @abc.abstractmethod
    def runComponent(self):
        raise NotImplementedError("Each class needs to implement runComponent")
        
    @abc.abstractmethod    
    def findFailedSamples(self):
        raise NotImplementedError("Each class needs to implement findFailedSamples")
    
#     def readLogFile(self, log_fn):
#         self.log = PlinkLogReader(log_fn)
        
    def getNoStartSamples(self):
        return self.log.getNoStartSamples()
        
    def getNoStartVariants(self):
        return self.log.getNoStartVariants()
    
    def getNoEndSamples(self):
        return self.log.getNoEndSamples()
        
    def getNoEndVariants(self):
        return self.log.getNoEndVariants()
        
    def getNoRemovedSamples(self):
        return self.log.getNoRemovedSamples()
        
    def getNoRemovedVariants(self):
        return self.log.getNoRemovedVariants()
    
    def writeFailedSamples(self, fail_fn = None):
        fn = self.failed_fn
        if fail_fn != None:
            fn = fail_fn
        with open(fn, 'w') as out_f:
            for fam_id, p_id in self.failed_samples:
                out_f.write("%s\t%s\n" % (fam_id, p_id))
        return fn
    
    def removeIndividuals(self, dataset_fn, failed_samples = None, output_fn = None):
        """
        Removes the persons listed in self.failed_samples file (can give other file name as option)
        Returns the name of the created file. 
        """
        fail_fn = self.failed_fn
        if failed_samples != None:
            fail_fn = failed_samples
        switch1 = "--remove %s" % fail_fn
        if output_fn == None:
            output_fn = "removed_%s" % (self.fail_type)
        if output_fn.find(self.args.output_dir) == -1:
            output_fn = os.path.join(self.args.output_dir, output_fn)
        self.plinkRunner.runPlinkCommand([switch1, "--make-bed"], dataset_fn, output_fn)
        self.tempFilesCreated(output_fn)
        return output_fn
        
    def tempFilesCreated(self, base_name):
        """
        Adds any files matching pattern basename.* to temp files (to be deleted later)
        """
        files = glob.glob(base_name + ".*")
        if len(files) == 0:
            print ("No files found for %s !!!" % base_name)
        for fn in files:
            self.temp_files.append(fn)
        
    def getTempFiles(self):
        return self.temp_files
    
    