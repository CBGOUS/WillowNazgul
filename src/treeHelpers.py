# -*- coding: utf-8 -*-
"""
Created on Fri Nov 19 14:58:51 2021

@author: norab


CONTAINS: 
    
    Help functions to ease data-loading.
    

"""


import pandas as pd
import os
from os.path import isfile, join
import re


def make_filelist(input_files):
    """
    Input: 
        input_files: path to folder
        type: str
        
    Function: 
        creates list of full path to every file in folder
        
    """
    try: 
        if isinstance(input_files, str):
            if '.csv' in input_files: 
                files = pd.read_csv(input_files, header = None)
                files = files[0].values.tolist()
            else: 
                files = [join(input_files, f) for f in os.listdir(input_files) 
                             if isfile(join(input_files, f))]
                files = [f.strip() for f in files]
            return files
    except Exception as e:
        print("Error with make_filelist:")
        print(e)

def geneName(phydists_file, filetype = "phydist", idtype = 'all'):
    """
    Input: 
        file: string filepath to distance matrix-file. 
        filetype: ninja tree or phydist file
        idtype: select gene id to return
    Function: 
        filters out name of gene the tree represents and assign to
        both Ensembl and gene name identifiers included on the form: 
            'ENSG00000000938___FGR'
    return: gene name (str)
            
    """
    
    if filetype == "phydist":
        subName = re.sub('^.*ENS', 'ENS', phydists_file)
        gene_name = re.sub('___CopD.csv$','', subName)
    else: 
        subName = re.sub('^.*ENS', 'ENS', phydists_file)
        gene_name = re.sub('_HUMAN.*','', subName)
    
    if idtype == 'ENS':
        gene_name = re.sub('___.*', '', gene_name)
    elif idtype == 'gene':
        gene_name = re.sub('^.*___','', gene_name)
    
    return gene_name

def load_phydists(file):
    """
    Input: folder: directory with file(s) 
    Function: read phydist file
    Return: distance matrix
    """
    return pd.read_csv(file, index_col = 0)

def load_simData(file='C:/Users/norab/Master/data/simulation/simNull_pvals.csv'):
    return pd.read_csv(file)

def load_uniqseqs(file= 'C:/Users/norab/Master/thesis_data/meta_data/9381_uniqseqs.txt'):
    """
    Input: file_path: str filedire containing number of unique sequences
        each tree is built from    
    Funtion: read input data
    Return: pandas dataframe-formatted input data
    """
    
    uniqseq = pd.read_csv(file, delimiter = ":", header = 0, index_col = None) 
    uniqseq.rename(index=lambda s: re.sub('_HUMAN__uniq.*', '', s), 
                   columns = {uniqseq.columns[0]: "uniqseq"}, inplace = True)
    uniqseq.reset_index(inplace=True)
    uniqseq.rename(columns = {"index": "gene"}, inplace=True)
    
    return uniqseq

def load_nz_phydists(file= 'C:/Users/norab/Master/thesis_data/meta_data/nonzero_phydists.csv'):
    """
    Input: file: str path to file with calculated non-zero phydists
    Funtion: Load data
    Return: Pd dataframe with values and gene names as indices.
    """
    
    totdist = pd.read_csv(file, "r", delimiter = ",", index_col = 0)
    totdist.rename(index=lambda s: re.sub('_HUMAN__full.*', '', s), 
                    columns = {totdist.columns[0]: "totdist"}, inplace = True)
    totdist.rename(index=lambda s: re.sub('^.*trees/', '', s), 
                    columns = {totdist.columns[0]: "totdist"}, inplace = True)
    totdist.reset_index(inplace=True)
    totdist.rename(columns={"index":"gene"}, inplace = True)
    
    return totdist

def load_GDRs(file = 'C:/Users/norab/Master/thesis_data/result_data/GDR/GDR_all.csv', group_category = 'all'):
    """ 
    Input: file: file with calculated GDR values
          group_category: option: 'super' and 'sub' 
    Return: pd DataFrame containing GDR values for super and/or 
    sub popultions for each gene
    
    """
    GDRs = pd.read_csv(file, header=0, index_col=False)
    
    if group_category == 'super':
        return GDRs[GDRs['level'] == 'super']
    elif group_category == 'sub':
        return GDRs[GDRs['level'] == 'sub']
    else:
        return GDRs 
    
def load_totdist_pops(file = 'C:/Users/norab/Master/data/other_measures/totdist_pops_all.csv'):
    """
    Input: file: file with total distance within each defined population
                for all genes.  
    Function: load data
    Return: pd dataframe with input data values 
    """

    return pd.read_csv(file)

def load_GDRnull(file_path = 'C:/Users/norab/Master/thesis_data/test_result_data/GDRnull_all_27.11.21.csv',
                 group_category = "All"):
    """ 
    Input: Filepath
    Return: pd DataFrame containing GDR values for super and sub popultions for each tree/gene
    
    """ 
    
    GDRs = pd.read_csv(file_path, header=0, index_col=False)
    
    if 'super' in group_category:
        return GDRs[GDRs['level'] == 'nullSuper']
    elif 'sub' in group_category:
        return GDRs[GDRs['level'] == 'nullSub']
    else:
        return GDRs 
        
    # YOU STOPPED HERE: WRANGLE so you can load this properly. All 
    
    
def get_gene_entry(gene, df):

    return df[df['gene'].str.contains(gene)]