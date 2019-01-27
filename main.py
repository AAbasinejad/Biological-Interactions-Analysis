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
    SGI_biogrid, SGI_iid = i.seed_gene_interactons(SG, biogrid, iid)
    union_biogrid, union_iid, all_biogrid_interactomes, all_iid_interactomes = i.union_interactions(SG, biogrid, iid)
    I_interactions = i.intersection_interactions(SG, biogrid, iid)
    SGI_integrated, union_integrated = i._integrations(SG, biogrid,iid)
    interactions_detail = i.detail_interactions_results(union_biogrid, union_iid, all_biogrid_interactomes, all_iid_interactomes)
    i.save_results(SGI_biogrid, SGI_iid, union_biogrid, union_iid, all_biogrid_interactomes, all_iid_interactomes, 
                                I_interactions, SGI_integrated, union_integrated, interactions_detail)
    
    # Netwrok analysis
    #union_adjacency, sgi_adjacency, I_adjacency = na.adjacency_matrix(SG, biogrid, iid) 
    #G_union, G_sgi, G_I, G_lcc_union, G_lcc_I = create_graph(union_adjacency, sgi_adjacency, I_adjacency)




if __name__ == '__main__':
    main()

