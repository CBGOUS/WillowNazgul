# -*- coding: utf-8 -*-
"""
Created on Fri May 21 08:32:27 2021

@author: norab
"""

import pandas as pd
import re
from random import shuffle


class treeInfo:
    
    """
    treeInfo class handles all input information and stores to class variables.
    
    """

    def __init__(self):
        
        """
        Function: 
            initiate class variables
        """
        self.name = None                # name of tree, eg. gene
        self.phydists = None            # phylogenetic distance
        self.categories = None          # categories for grouping structures
        self.group_info = None          # groups in each category
        self.sample_info = None         # info of each sample in tree
        self.random_groups = False      # if samples are randomly assigned to groups
        
    def setup(self, phydists_file=None, group_info_file=None, categories=None):
        """
        Input: 
            phydists_file: 
                filepath to distance matrix (CSV)
                type: string
            
            group_info_file: 
                filepath group category information (CSV)
                type: string
                
        Function:
            call other treeInfo-functions to format and store input information
        """
        
        if phydists_file:
            self.setPhydists(phydists_file)
            self.setName(phydists_file)
            
        if categories != 'None': 
            self.setCategories(categories)
        
        if group_info_file != 'None':
            self.setGroupInfo(group_info_file)
        
        
        if (phydists_file != 'None') and (categories != 'None'):
            self.setSampleInfo()

    def setPhydists(self, phydists_file):
        """
        Input:
            phydists_file: filepath to distance matrix (CSV)
                type: str
        
        Function: 
            read distance matrix values, assign information to class variable
        """
        self.phydists = pd.read_csv(phydists_file, index_col = 0, 
                                    dtype={'a': str})
    
    def setName(self, phydists_file):
        """
        Input: 
            phydists_file: filepath to distance matrix (CSV)
                type: str
        
        Function: 
            read name of distance matrix file, assign to class variable
            
        Note: 
            hard-coded for customized format used in master project.            
            both ENSEBL and gene name identifiers included on the form
            'ENSEBL___GENENAME'
            example: 'ENSG00000000938___FGR'
            
        TO UPDATE
            
        """

        subName = re.sub('^.*ENS', 'ENS', phydists_file)
        self.name = re.sub('___CopD.csv$','', subName)
            
    def setCategories(self, categories, sep = '___'):
        """
        Input: 
            categories: group category information
            type: str
                format: 'category1___category2___category3___...____'
                example: 'SUPER___SUB'
            
            sep: separator of category names in string, default: ___
            
        Function: 
            assign categories class variable

        """
        self.categories = categories.split(sep)
    
    def setSampleInfo(self):
        """ 
        Function:
            - assign sample information to class variable
            - samples are stored in nested dictionary
            - keys to nested dictionary indices are equal to numeric order of
            samples in rows and columns of distance matrix
                     
            dict format example: 
                {0: {'sample': 'HG02545', 'SUPER': 'AFR', 'SUB': 'ACB'}
                1: {'sample': 'HG00334', 'SUPER': 'EUR', 'SUB': 'FIN'}}
        
        """
        self.sample_info = {}
        
        for ind, sample_in in enumerate(self.phydists.columns):
            sample_in = sample_in.split('___')
            # Add sample info
            sample_trans = {'sample':sample_in[-1]}
            
            # Add all categories
            for cat_ind, cat in enumerate(self.categories):
                sample_trans[cat] = sample_in[cat_ind]
    
            # Store sample info
            self.sample_info[ind] = sample_trans

    def shuffleSampleInfo(self):
        """
        Function: 
            shuffle sample information in each category. rows and column 
            indices of distance matrix will no longer correspond to original 
            sample informaiton

        """

        for cat in self.categories:
            samples = [sample[cat] for sample in self.sample_info.values()]
            shuffle(samples)
            
            for ind in self.sample_info:
                new = self.sample_info[ind]
                new[cat] = samples.pop()
                self.sample_info[ind] = new

        self.random_groups = True

    def setUniqseqMap(self, uniqseq_map_file):
        """
        Input: 
            uniqseq_map_file: filepath containing unique sequences
            type: str
            
        Function: 
            retrieve unique sequence informaiton
        
        """
        # Read file
        file = pd.read_csv(uniqseq_map_file, header = None)
        
        # Counts
        self.uniqseq_count = len(file[0])   # What is this for
        self.uniqseq_map = sorted(file.count(axis=1), reverse=True)

    def setGroupInfo(self, group_info_file):
        """
        Input: 
            group_info_file: filepath to file with group info (CSV)
            type: str
            
        Function: 
            organize information into a list containig a set of defined 
            groups for each group category 
            example: 
                two categories: "super" and "sub"
                list = [{AFR, EUR, EAS}, {FIN, YRI, GBR}]
        """
        group_info = pd.read_csv(group_info_file, delimiter='\t')
        
        self.group_info = dict()
        for cat in self.categories:
            groups = group_info[group_info['GroupCategory'] == cat]
            groups = set(groups['GroupName'])
            self.group_info[cat] = groups

    def getGroupInfo(self): return self.group_info
    def getSampleInfo(self): return self.sample_info
    def getName(self): return self.name
    def getPhydists(self):return self.phydists
    def getUniqseqMap(self): return self.uniqseq_map
    def getUniqseqCount(self): return self.uniqseq_count
    def getCategories(self): return self.categories
    

if __name__ == '__main__':
  
    # Test with simple case
    phydists_file = 'C:/Users/norab/Master/WillowProject/Willow1.0/jobs/testOneGene/job_input/geneDists/ENSG00000000938___FGR___CopD.csv'
    group_info_file = 'C:/Users/norab/Master/WillowProject/Willow1.0/jobs/testOneGene/job_input/phydist_population_classes.tsv'
    categories = 'SUPER___SUB'
    test_tree = treeInfo()
    test_tree.setup(phydists_file, group_info_file, categories)
    
    sample_info = test_tree.getSampleInfo()
    group_info = test_tree.getGroupInfo()
    categories = test_tree.getCategories()


