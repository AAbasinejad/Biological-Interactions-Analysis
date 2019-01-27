#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import Basic_info as Bi
import Network_Analysis as na
import Interactions as i
import pandas as pd
import sys



def print_usage():
    print('usage: python main.py seed_gene.txt biogrid_dataset iid_dataset')
    print('seed_gene file: The file that contains only seed gene symbols')
    print('biogrid_dataset: The latest HOMO Sapiens BIOGRID data delimited in delimited by tab')
    print('iid_dataset: The latest iid_human dataset')

def check_input_style(input_list):
    try:
        SG, biogrid, iid = read_input(input_list)
    # if no input is given, print out a usage message and exit
    except:
        print_usage()
        sys.exit(0)
        return 
    return SG, biogrid, iid
def read_input(input_list):
    if (len(input_list)!= 4):
        return 0 
    else:
        seed_gene_file = open(input_list[1], 'r')
        seed_genes = seed_gene_file.read().split(',')
        SG = seed_genes
        seed_gene_file.close()
        biogrid = pd.read_csv(input_list[2], delimiter = '\t', 
                              usecols = ['Official Symbol Interactor A', 'Official Symbol Interactor B'] ).assign(uniprot1 = 'NA', uniprot2 = 'NA', dbs = 'BioGrid')
        biogrid = biogrid.rename(columns={'Official Symbol Interactor A':'interactor A',
                                          'Official Symbol Interactor B':'interactor B'})
        iid = pd.read_csv(input_list[3], delimiter=' |\t', 
                          usecols=['symbol1', 'symbol2', 'uniprot1', 'uniprot2', 'dbs'])
        iid = iid.rename(columns={'symbol1':'interactor A', 'symbol2': 'interactor B'})

    
    return SG, biogrid, iid
    



def main():    
    
    input_list = sys.argv
    SG, biogrid, iid = check_input_style(input_list)

    # Basic Info

    Bi.Basic_informations(SG, iid)

    # Interactions

    


    # Netwrok analysis
    union_adjacency, sgi_adjacency, I_adjacency = na.adjacency_matrix(SG, biogrid, iid) 





if __name__ == '__main__':
    main()

