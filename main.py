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

def read_input(input_list):
    if (len(input_list)!= 4):
        print_usage()
        return
    else:
        seed_gene_file = open(input_list[1], 'r+')
        SG = seed_gene_file.read().split('\n')
        seed_gene_file.close()
        biogrid = pd.read_csv(input_list[2], delimiter = '\t', 
                              usecols = ['Official Symbol Interactor A', 'Official Symbol Interactor B'] ).assign(uniprot1 = 'NA', uniprot2 = 'NA', dbs = 'BioGrid')
        biogrid = biogrid.rename(columns={'Official Symbol Interactor A':'interactor A',
                                          'Official Symbol Interactor B':'interactor B'})
        iid = pd.read_csv(input_list[3], delimiter=' |\t', 
                          usecols=['symbol1', 'symbol2', 'uniprot1', 'uniprot2', 'dbs'])
        iid = iid.rename(columns={'symbol1':'interactor A', 'symbol2': 'interactor B'})
        iid = iid[['interactor A', 'interactor B', 'uniprot1', 'uniprot2', 'dbs']]
    
    return SG, biogrid, iid
    
#['main.py', 'seed_genes.txt', 'BIOGRID-ORGANISM-Homo_sapiens-3.5.166.tab2.txt', 'iid.human.2018-05.txt']
def main():    
    
    input_list = sys.argv
    SG, biogrid, iid = read_input(input_list)
    

    # Basic Info

    #Bi.Basic_informations(SG, iid, biogrid)

    # Interactions
    SGI_biogrid, SGI_iid = i.seed_gene_interactons(SG, biogrid, iid)
    union_biogrid, union_iid, all_biogrid_interactomes, all_iid_interactomes = i.union_interactions(SG, biogrid, iid)
    I_interactions = i.intersection_interactions(SG, biogrid, iid)
    SGI_integrated, union_integrated = i._integrations(SG, biogrid,iid)
    interactions_detail = i.detail_interactions_results(union_biogrid, union_iid, all_biogrid_interactomes, all_iid_interactomes, SG)
    i.save_results(SGI_biogrid, SGI_iid, union_biogrid, union_iid, all_biogrid_interactomes, all_iid_interactomes, 
                   I_interactions, SGI_integrated, union_integrated, interactions_detail, SG)
    
    # Netwrok analysis
    #union_adjacency, sgi_adjacency, I_adjacency = na.adjacency_matrix(SG, biogrid, iid) 
    #G_union, G_sgi, G_I, G_lcc_union, G_lcc_I = create_graph(union_adjacency, sgi_adjacency, I_adjacency)




if __name__ == '__main__':
    main()

