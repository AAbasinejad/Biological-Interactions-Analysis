#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import os




# all interactions that involves seed genes only
def seed_gene_interactons(SG, biogrid, iid):
    SGI_biogrid = biogrid.loc[biogrid['interactor A'].isin(SG) & 
                              biogrid['interactor B'].isin(SG)]
    SGI_biogrid = SGI_biogrid.drop_duplicates(subset=['interactor A','interactor B'])
    
    SGI_iid = iid.loc[iid['interactor A'].isin(SG) & 
                      iid['interactor B'].isin(SG)]
    SGI_iid = SGI_iid.drop_duplicates(subset=['interactor A', 'interactor B'])
   
    return SGI_biogrid, SGI_iid

# union_interaction is an interaction that at least one seed_gene involves  
def union_interactions(SG, biogrid, iid):
    # interactions which which at least one seed_gene involves
    biogrid_sg_interactions = biogrid.loc[biogrid['interactor A'].isin(SG) | 
                                          biogrid['interactor B'].isin(SG)]
    all_directed_biogrid_interactomes = set(biogrid_sg_interactions["interactor A"]).union(set(biogrid_sg_interactions["interactor B"]))
    # all interactions with consideration of interactomes interactions
    union_biogrid = biogrid.loc[biogrid['interactor A'].isin(all_directed_biogrid_interactomes) &
                                biogrid['interactor B'].isin(all_directed_biogrid_interactomes)]
    union_biogrid = union_biogrid.drop_duplicates(subset=['interactor A','interactor B'])
    
    #####
    iid_sg_interactions = iid.loc[iid['interactor A'].isin(SG) |
                                  iid['interactor B'].isin(SG)]
    all_directed_iid_interactomes = set(iid_sg_interactions['interactor A']).union(set(iid_sg_interactions['interactor B']))
    #####
    union_iid = iid.loc[iid['interactor A'].isin(all_directed_iid_interactomes) &
                        iid['interactor B'].isin(all_directed_iid_interactomes)]
    union_iid = union_iid.drop_duplicates(subset=['interactor A','interactor B'])

    return union_biogrid, union_iid, all_directed_biogrid_interactomes, all_directed_iid_interactomes

# all proteins interacting with at least one seed gene confirmed by both DBs 
def intersection_interactions(SG, biogrid, iid):
    u_bio, u_iid, all_bio, all_iid = union_interactions(SG, biogrid, iid)
    
    I_interactions = pd.merge(u_iid, u_bio, how='inner', on=['interactor A', 'interactor B'])
    
    return I_interactions
# Integrating the results of SGI and biogrid parts, (the interested result is this one)
def _integrations(SG, biogrid, iid):
    SGI_biogrid, SGI_iid = seed_gene_interactons(SG, biogrid, iid)
    u_bio, u_iid, all_bio, all_iid = union_interactions(SG, biogrid, iid)
    
    SGI_frames = [SGI_biogrid, SGI_iid]
    SGI_integrated = pd.concat(SGI_frames)
    SGI_integrated = SGI_integrated.drop_duplicates(subset=['interactor A', 'interactor B'])
    
    
    union_frames = [u_bio, u_iid]
    union_integrated = pd.concat(union_frames)
    union_integrated = union_integrated.drop_duplicates(subset=['interactor A', 'interactor B'])
    
    return SGI_integrated, union_integrated

# making a table to show some stats about two databases - seperatedly
def detail_interactions_results(u_bio, u_iid, all_bio, all_iid, SG):
    
    db_names = ['BioGrid', 'IID']
    Num_of_founded_SG = list() 
    total_interacting_proteins = list()
    total_interactions = list()
    
    Num_of_seed_genes_biogrid = (set(u_bio['interactor A']).union(set(u_bio['interactor B']))).intersection(set(SG))
    Num_of_seed_genes_iid = (set(u_iid['interactor A']).union(set(u_iid['interactor B']))).intersection(set(SG)) 
    Num_of_founded_SG = [len(Num_of_seed_genes_biogrid), len(Num_of_seed_genes_iid)]
    
    total_interacting_proteins = [len(all_bio), len(all_iid)]
    total_interactions = [u_bio.shape[0], u_iid.shape[0]]
    d = {'DB_name': db_names, 'Num. of founded SG': Num_of_founded_SG, 
         'Total interacting Prot.': total_interacting_proteins, 'total interactions': total_interactions}
    interactions_detail = pd.DataFrame(data=d)
    
    return interactions_detail

# saving results in csv and txt formats...
def save_results(SGI_bio, SGI_iid, u_bio, u_iid, all_bio, all_iid, I_interactions, SGI_integrated, union_integrated, interactions_detail, SG):
        
    os.mkdir('Interactions_Results', mode = 0o777)
    os.chdir('Interactions_Results')
    
    # save the results in csv files 0 SGI
    SGI_bio.to_csv('SGI_biogrid.csv', sep='\t')
    SGI_iid.to_csv('SGI_iid.csv', sep='\t')
    # write details of results in a text file
    f = open('interactions_detailed_result.txt', 'w')
    f.write('seed_gene_interactions(a.k.a SGI) results....\n')
    f.write('BioGrid data set SGI Interactions: '+ str(SGI_bio.shape[0]) + '\n')
    f.write('iid data set SGI Interactions: '+ str(SGI_iid.shape[0]) + '\n')
    
    # save the results in csv files - union
    u_bio.to_csv('Union_biogrid.csv', sep='\t')
    u_iid.to_csv('Union_iid.csv', sep='\t')
    # write details of results in a text file
    f.write('union_interactions results....\n')
    f.write('BioGrid data set directed interactoms: '+ str(len(all_bio - set(SG))) + '\n')
    f.write('BioGrid data set union_interactions: '+ str(u_bio.shape[0]) + '\n')
    f.write('iid data set directed interactoms: '+ str(len(all_iid - set(SG))) + '\n')
    f.write('iid data set union_interactions: '+ str(u_iid.shape[0]) + '\n')
    
    # intersection results
    I_interactions.to_csv('Intersection_interaction.csv', sep='\t')
    f.write('Intersection_interactions results....\n')
    f.write('Number of interactions confirmed by both DBs: '+ str(I_interactions.shape[0])+ '\n')
    
    # Integration results
    SGI_integrated.to_csv('SGI_integrated.csv', sep='\t')
    union_integrated.to_csv('union_integrated.csv', sep='\t')
    f.write('integration results....\n')
    f.write('total SGI interactions: '+ str(SGI_integrated.shape[0])+ '\n')
    f.write('total union interactions: '+ str(union_integrated.shape[0])+ '\n')
    
    # interactions_details
    interactions_detail.to_csv('interactions_detail_results.csv', sep='\t')
    
    f.close()
    
    return
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    