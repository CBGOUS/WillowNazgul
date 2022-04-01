# -*- coding: utf-8 -*-
"""
Created on Thu May  6 11:33:03 2021

@author: norab

"""

import pandas as pd
from treeInformation import treeInfo

class treeMetrics(treeInfo):
    
    """
    inherits treeInfo class. 
    class to perform calculations utilizing informaiton stored with treeInfo. 
    """
    
    def __init__(self):
        """
        Function: 
            - initiate class variables
            - inherit treeInfo
        """
        
        self.group_dists = None
        self.mean_group_dists = None
        self.GDRs = None
        self.nz_phydists = None
                
        treeInfo.__init__(self)
        
    def calcGroupDists(self):
        """
        Function: 
            calculate inter- and intra group distances for all categories
        """

        phydists = self.phydists.to_numpy()

        catWithSums = dict()
        catBetSums = dict()
        
        for category, groups in self.group_info.items():
            groupWithSums = dict()
            groupBetSums = dict()
            
            for group in groups: 
                groupWithSums[group] = [0,0]
                groupBetSums[group] = [0,0]
                
            catWithSums[category] = groupWithSums
            catBetSums[category] = groupBetSums
        
        row_length = len(phydists)
        row_start = 1
        
        for sample1 in range(row_length):
            
            for sample2 in range(row_start, row_length):
                dist_val = phydists[sample2][sample1]        # Distance value
                
                for cat in self.categories: 
                    g1 = self.sample_info[sample1][cat]
                    g2 = self.sample_info[sample2][cat]
                    
                    if g1 == g2:
                        catWithSums[cat][g1] = [old+new for old, new in zip(catWithSums[cat][g1], [dist_val,1])]
                    else: 
                        catBetSums[cat][g1] = [old+new for old, new in zip(catBetSums[cat][g1], [dist_val,1])]
                        catBetSums[cat][g2] = [old+new for old, new in zip(catBetSums[cat][g2], [dist_val,1])]
            
            row_start += 1  
                
        self.group_dists = {'catWithSums': catWithSums, 'catBetSums': catBetSums}

    def calcMeanGroupDists(self):
        """
        Function: 
            calculates mean distance between all pairwise samples within- and
            between groups respectivly, for each category
        """
        
        if not self.group_dists: self.calcGroupDists()
        
        catWithSums = self.group_dists['catWithSums']
        catBetSums = self.group_dists['catBetSums']
        
        self.mean_group_dists = {}
        for cat in self.categories:
            
            withSums = catWithSums[cat]
            betSums = catBetSums[cat]
            
            dist_summary = {'with':0, 'bet':0}
            count_summary = {'with':0, 'bet':0}
            
            for group, val in withSums.items(): 
                dist_summary['with'] += val[0]
                count_summary['with'] += val[1]
            
            for group, val in betSums.items(): 
                dist_summary['bet'] += val[0]
                count_summary['bet'] += val[1]
            
            mean_dists = {key: round(dist_summary[key] / count_summary[key], 6)  # Was 8 
                          for key, val in dist_summary.items() if count_summary[key]}
          
            self.mean_group_dists[cat] = mean_dists
        
    def calcGDR(self):
        """
        Function: 
            calculate GDR
        """
        
        if not self.mean_group_dists: self.calcMeanGroupDists()
            
        self.GDRs = {}
        for cat, val in self.mean_group_dists.items():
            
            if val['bet']:
                GDR = round(val['with'] / val['bet'], 4)
            else: 
                GDR = 1     # Maybe should be NAN - think
            
            self.GDRs[cat] = GDR
            
    def calcNonZeroPhydists(self):
        """
        Input:
            percent: give nonero variable as percentage
    
        Function:
            calculate proportion of distance matrix that has non-zero values, 
            assign to class variable
        """

        nonZero_row = pd.DataFrame((self.phydists != 0).astype(int).sum(axis=1))
        nonZeroCount = int(nonZero_row.sum(axis=0))
        
        num_entries = (self.phydists.shape[0] * self.phydists.shape[1]) - self.phydists.shape[0]
        
        self.nz_phydists = nonZeroCount / num_entries
        
    def getGroupDists(self): 
        if not self.group_dists: self.calcGroupDists()
        return self.group_dists
    
    def getMeanGroupDists(self):
        if not self.mean_group_dists: self.calcMeanGroupDists()
        return self.mean_group_dists
        
    def getGDRs(self): 
        if not self.GDRs: self.calcGDR()
        return self.GDRs
    
    def getNZphydists(self):
        if not self.nz_phydists: self.calcNonZeroPhydists()
        return self.nz_phydists


if __name__ == '__main__':
    pass
    # Test with simple case
    
    # # Test files
    # phydists_file = 'C:/Users/norab/Master/WillowProject/Willow1.0/jobs/testOneGene/job_input/geneDists/ENSG00000000938___FGR___CopD.csv'
    # group_info_file = 'C:/Users/norab/Master/WillowProject/Willow1.0/jobs/testOneGene/job_input/phydist_population_classes.tsv'
    # categories = 'SUPER___SUB'
    
    
    # test_tree = treeMetrics()
    # test_tree.setup(phydists_file, group_info_file, categories)
    
    # sample_info = test_tree.getSampleInfo()
    # group_info = test_tree.getGroupInfo()
    # categories = test_tree.getCategories()
    
    # # Calc group dists
    # starttime = timeit.default_timer()
    # timeGroupDists = test_tree.calcGroupDists()
    # endtime = timeit.default_timer()
    
    # test_tree.calcMeanGroupDists()
    # test_tree.calcGDR()    
    
    # GDRss = test_tree.getGDRs()