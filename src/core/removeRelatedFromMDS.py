'''
Created on 18.8.2016

@author: paakkone
'''

import os
import random
from component import Component

class RemoveRelatedFromMDS(Component):
    
    def __init__(self, args, failed_fn = None, plink_commands_f = None):
        self.fail_type = "mds_rel"
        super().__init__(args, failed_fn = failed_fn, plink_commands_f = plink_commands_f)
        
    def addIDCount(self, id_counts, fid, iid):
        new_id = fid + "&" + iid
        if new_id not in id_counts.keys():
            id_counts[new_id] = 0
        id_counts[new_id] += 1
        return id_counts
    
    def getIDCounts(self, pihat_rows):
        """
        Counts the occurances of the ID's
        """
        id_counts = {}
        for fid1, iid1, fid2, iid2, pi_hat in pihat_rows:
            id_counts = self.addIDCount(id_counts, fid1, iid1)
            id_counts = self.addIDCount(id_counts, fid2, iid2)
        return id_counts
    
    def removeSelectedIDRows(self, pihat_rows, selected_id):
        """
        Delete rows from pihat_rows which contains the selected ID
        """
        remove_fid, remove_iid = selected_id.split("&")
        # don't actually delete items but create a new list with those rows which don't have the ID (simpler than deleting)
        remaining_rows = []
        #for i, fid1, iid1, fid2, iid2, pi_hat in enumerate(pihat_rows):
        for i, row in enumerate(pihat_rows):
            fid1, iid1, fid2, iid2, pi_hat = row
            if not (fid1 == remove_fid and iid1 == remove_iid) and not (fid2 == remove_fid and iid2 == remove_iid):
                remaining_rows.append(pihat_rows[i])
        return remaining_rows    
    
    def runPlink(self):
        # plink not needed at this step
        pass
    
    def runComponent(self, dataset_fn, pi_val_fn):
        # filter for failures
        self.findFailedSamples(pi_val_fn, dataset_fn)
        # write out failures
        self.writeFailedSamples()
        # remove samples and write new file set
        new_dataset_fn = self.removeIndividuals(dataset_fn)
        # read log file for reporting of number of samples
        self.log.readLogFile(new_dataset_fn + ".log")
        # return new file set name
        return new_dataset_fn
    
    
    def findFailedSamples(self, pi_val_fn, dataset_fn):
        """
        Reads in the pihat value file which contains all sample pairs with pihat values greater than given limit.
        For MDS these samples need to be removed (or half of them). For that it checks for repeated ID's and
        removes those which are several times in the file first, until each sample is just once in the file.
        Then from those one sample per row is randomly selected for removal.
        """
        delim = "&"
        fail_type = "mds_rel"
        pihat_rows = []
        # id list of id's which will be removed
        ids_to_remove = []
        header = True
        with open(pi_val_fn) as in_f:
            for line in in_f:
                if header:
                    header = False
                else:
                    words = line.split()
                    pihat_rows.append(words)
                    if len(words) != 5:
                        print ("*********Only four entries! %s" % words)
        # values read, now go through rows and remove IDs from it starting with the most common id's
        done = False
        if len(pihat_rows) == 0:
            # skip this step if there were no pihat values over the treshold in the data
            done = True
        while not done:
            id_counts = self.getIDCounts(pihat_rows)
            max_repetitions = max(id_counts.values())
            if max_repetitions > 1:
                # find all max value ids
                max_val_sample_ids = [sid for sid in id_counts if id_counts[sid] == max_repetitions]
                # choose randomly one of the max repetition value ID's
                selected_id = random.choice(max_val_sample_ids)
            else:
                # go through rest of the lines and choose from each line randomly one key
                fid1, iid1, fid2, iid2, pi_hat = pihat_rows[0]
                selected_id = fid1 + delim + iid1
                if random.randint(1,2) == 2:
                    selected_id = fid2 + delim + iid2                    
            # put selected value to ids_to_remove list
            ids_to_remove.append(selected_id)
            # remove rows which contain selected id
            pihat_rows = self.removeSelectedIDRows(pihat_rows, selected_id)
            # check if there is any rows left
            if len(pihat_rows) == 0:
                done = True
        # write out to_remove list
        for sid in ids_to_remove:
            fid, iid = sid.split(delim)
            self.failed_samples.append([fid, iid])
        return
