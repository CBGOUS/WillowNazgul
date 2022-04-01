# -*- coding: utf-8 -*-
"""
Created on Wed Jun  2 05:48:23 2021

@author: norab
"""

import csv
from datetime import datetime
import pandas as pd
from treeMetrics import treeMetrics
from treeHelpers import make_filelist
import logging

class treeRun:
    """
        class to organize tree instance calculations for input files
    """
    
    def __init__(self, config):
        
        """
        Function: 
            initiate class variables
        """
        
        # Initialize logger
        logger = logging.getLogger(__name__)
        logger.info("--initialising")
        
        # setup configuration
        self.config = config
        self.config_setup()
		# create list of phylogenetic distance files to run
        self.file_list = make_filelist(input_files=self.input_phydist_folder)
        
        # Run! 
        self.run_func()

    def config_setup(self):
        
        """
        Input: 
            config: configuration filepath  (yaml)
            type: str
                       
        Function: 
            configure variables from file
        
        """
        
        # Make datetime stamp
        for key, val in self.config.items():
            if 'datetime' in str(val): 
                self.config[key] = val.replace('datetime', datetime.now().strftime("%d.%m.%Y_%H.%M"))
        
        # Folder/file paths
        #self.input_folder = config.get('input_folder').strip()
        self.input_phydist_folder = self.config.get('input_phydist_folder').strip()
        self.output_folder = self.config.get('output_folder').strip()
        self.output_unprocessed = self.output_folder + self.config.get('output_unprocessed').strip()
               
        # Input information
        self.func = self.config.get('func').strip()
        self.group_info = self.config.get('input_group_info').strip()
        self.filter_select = self.config.get('filter_select').strip()
        self.filter_skip = self.config.get('filter_skip').strip()
        self.group_categories = self.config.get('group_categories').strip()
        self.num_random_values = int(self.config.get('num_random_values'))
        self.random = self.num_random_values

    def run_calcGDR(self):
        """
        Funcition: 
            - calculate GDR value for all given input distance matrix files
            - save values to file
            - if error with file, filepath is saved to separate file and
            calculation proceeds with next file
        """
        
        try:        
            if '.csv' in self.filter_select:
                filter_select = pd.read_csv(self.filter_select, header=None)
                filter_select.columns = ['name']
                filter_select= list(filter_select['name'])
                
                keep_files = []
                for file in self.file_list:
                    if any(gene in file for gene in filter_select):
                        keep_files.append(file)
                        self.file_list = keep_files
            
            if '.csv' in self.filter_skip:
                filter_skip= pd.read_csv(self.filter_skip, header=None)
                filter_skip.columns = ['name']
                filter_skip= list(filter_skip['name'])
                
                remove_files = []
                for file in self.file_list:
                    if any(gene in file for gene in filter_skip):
                        remove_files.append(file)
        
                self.file_list = [f for f in self.file_list if f not in remove_files]
            
            ind = 0
            ind_len = len(self.file_list)
            
            if self.random: 
                num_iter = int(self.num_random_values)
            else:
                num_iter = 1

            for i in range(num_iter):
                
                for dist_file in self.file_list:
                    
                    
                    try:         
                        ind +=1
                        print(f"Processing file {ind} / {ind_len} ")
    
    
                        tree = treeMetrics()
                        tree.setup(dist_file.strip(), self.group_info, self.group_categories)
                        
                        if self.random: 
                            tree.shuffleSampleInfo()
    
                        tree.calcGDR()
                        GDRs = tree.getGDRs()
                        
                        for cat in GDRs:
                            
                            GDR_save = [tree.getName(), GDRs[cat]]
                            if self.random: 
                                save_to = self.output_folder + 'GDR_random_' + cat + '_' + str(datetime.now().strftime("%d.%m.%Y")) + '.csv'
                            else:
                                save_to = self.output_folder + 'GDR_' + cat + '_' + str(datetime.now().strftime("%d.%m.%Y")) + '.csv'
                        
                        
                            with open(save_to, 'a', newline='') as f:   # write to file    
                                writer = csv.writer(f)
                                writer.writerow(GDR_save)

        
                    except Exception as e: 
                        print("Error to prosess dist file:", dist_file)
                        print("Error: ", e)
                        # Save filepath to unprocessed file
                        file = str(dist_file)
                        with open(self.output_unprocessed, 'a') as f: 
                            f.write(file)

                        
                        # continue with next
                        pass
                
        except Exception as e: 
            print("error with runCalc:\n")
            print(e)
    
    def run_calcNZphydists(self):
        """
        Abb:      calculate non-zero phydists
        Input:    folder with phylgenetic distance files
        Function: calculate total amount of non-zero values for each file in folder
        Return:   pandas dataframe with values for all genes
        
        """
        
        # initiate df 
        all_nz_phydists = pd.DataFrame(columns = ['gene', 'nz_phydists'])
        
        for phydist_file in self.file_list:
            ind = 0
            ind_len = len(self.file_list)
            
            try:         
                ind +=1
                print(f"Processing file {ind} / {ind_len} ")
                
                tree = treeMetrics()
                tree.setup(phydist_file.strip(), self.group_info, self.group_categories)
                nz_phydists = tree.getNZphydists()
                geneID = tree.getName()
                
                all_nz_phydists = all_nz_phydists.append({'gene': geneID, 'nz_phydists':nz_phydists}, ignore_index = True)
                
                ind += 1
                
            except Exception as e: 
                print("Error to calc NZ phydist for file:" , phydist_file)
                print(e)
        
        

        save_to = self.output_folder + 'NZphydists' + '_' + str(datetime.now().strftime("%d.%m.%Y")) + '.csv'
           
        with open(save_to, 'a', newline='') as f:   # write to file    
            writer = csv.writer(f)
            writer.writerow(all_nz_phydists)

    def run_calcTest(self):
        """
        Test function to check functionality of program.
        """
        return 

    def run_func(self):
        """
        Run specified calculation procedure

        """
        if self.func == "calcGDR":
            self.run_calcGDR()
        elif self.func == "calcGDRrandom":
            self.random = True
            self.run_calcGDR()
        elif self.func == "calcNZphydists":
            self.run_calcNZphydists()
        elif self.func == "calcTest":
            return self.run_calcTest()
        else:
            print("No valid function called ", self.func)
   
        #logger.info('Session finished successfully: {0}'.format(datetime.now().strftime("%d.%m.%Y_%H.%M")))


import timeit

# if __name__ == '__main__':
#     pass


