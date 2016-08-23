'''
Created on 18.8.2016

@author: paakkone
'''
import sys
import os
import subprocess
import numpy

class PlinkRunner(object):
    def __init__(self, plink_commands_f = None):
        self.plink_path = "/apps/genetics/bin/plink-1.9"
        self.plink_commands_f = plink_commands_f

    def runPlinkCommand(self, switches, in_fn, out_fn):
        """
        builds and runs plink command
        """
        cmd = self.buildPlinkCommand(switches, in_fn, out_fn)
        if cmd == None:
            print ("Could not run plink, exiting program")
            sys.exit()
        self.runCommand(cmd)
        return True    
            
    def buildPlinkCommand(self, switches, in_fn, out_fn):
        """
        Gets switches for plink command and returns the command string
        """
        for ending in [".bim", ".bed", ".fam"]:
            fn = in_fn + ending
            if not os.path.exists(fn):
                print ("Plink input file %s does not exist!" % fn)
                sys.exit()
            if not os.path.isfile(fn):
                print ("Plink input file %s is not a file!" % fn)
                sys.exit()
        # start command and input file
        cmd = "%s --bfile %s" % (self.plink_path, in_fn)
        for sw in switches:
            cmd += " %s" % sw
        # output file
        cmd += " --out %s" % out_fn
        return cmd
        
    def runCommand(self, cmd):
        if os.name == "posix":
        # run SNPTest only in linux, otherwise unit tests will always fail
            print ("Running plink: %s" % cmd)
            if self.plink_commands_f != None:
                self.plink_commands_f.write("%s\n" % cmd)
            p = subprocess.Popen(cmd.split(), stdin = subprocess.PIPE, stdout = subprocess.PIPE, close_fds = True)
            p.wait()    
            
def getUpperAndLowerLimits(vals, std_limit):
    """
    Calculates upper and lower limit when given values and standard deviation limit
    """
    mean = numpy.mean(vals)
    std = numpy.std(vals, ddof=1)
    lower = mean - (std_limit * std)
    upper = mean + (std_limit * std)
    return lower, upper
