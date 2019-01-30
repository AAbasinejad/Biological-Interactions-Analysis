#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import pandas as pd
import requests
from bs4 import BeautifulSoup as bs
import warnings
warnings.filterwarnings("ignore")
import os


# obtaining Basic informations from web and local datasets
def Basic_informations(SG, iid, biogrid):
    os.mkdir('basic_info', mode = 0o777)
    os.chdir('basic_info')
    gene_ids = list()
    uniprot_AC = list()
    official_gene_symbol = list()
    official_protein_name = list()
    description = list()
    print('Fetching Basic Informations...')
    # the first loop is to get the gene IDs from ncbi.nlm.nih.gov which is also approved by HGNC
    for i in range(len(SG)):
        print(i)
        ids = list()
        em_list = list()
        index_list = list()
        response = requests.get('https://www.ncbi.nlm.nih.gov/gene/?term=' + SG[i])
        html = response.content
        soup = bs(html, 'html.parser')
        for items in soup.findAll('em'):
            em_list.append(items.text[0:4])
            for j in range(len(em_list)):
                if (em_list[j] == 'Homo'):
                    index_list.append(j)
        for items in soup.findAll('span', {'class', 'gene-id'}):
            ids.append(items.text[4:])
        gene_ids.append(ids[index_list[0]])
        
    # second loop is to get the proper Uniprot_AC of each seed gene by processing local iid data set 
    not_iid = open('Not_existed_seed_genes_iid.txt', 'w')
    for i in range(len(SG)):
        if (iid[iid['interactor A'] == SG[i]].shape[0] > 0):
            uniprot_AC.append(iid[iid['interactor A'] == SG[i]].iloc[0][2])
        elif (iid[iid['interactor B'] == SG[i]].shape[0] > 0):
            uniprot_AC.append(iid[iid['interactor B'] == SG[i]].iloc[0][3])
        else:
            uniprot_AC.append(SG[i])
            # seed_genes that does not exist in iid data_set
            not_iid.write(SG[i]+ '\n')
    not_iid.close()
    not_bio = open('Not_existed_seed_genes_biogrid.txt', 'w')
    for i in range(len(SG)):
        if (biogrid[biogrid['interactor A'] == SG[i]].shape[0] == 0 | 
            biogrid[biogrid['interactor B'] == SG[i]].shape[0] == 0):
            not_bio.write(SG[i]+ '\n')
    not_bio.close()
    # third (and final) loop is to go through all the fetched gene_ids in first part (open the proper webpage on ncbi.nlm.nih.gov) and fetch the interested
    # informations: official_gene_symbol, official_protein_name, breif_description
    for i in range(len(gene_ids)):
        print(i)
        desc_list = list()
        response = requests.get('https://www.ncbi.nlm.nih.gov/gene/'+ gene_ids[i])
        html = response.content
        soup = bs(html, 'html.parser')
        for items in soup.findAll(['dd']):
            desc_list.append(items.text)
        official_gene_symbol.append(desc_list[0][0:-16])
        official_protein_name.append(desc_list[1][0:-16])
        description.append(desc_list[9].split('.')[0]) 

    d = {'seed_genes': SG,'official gene symbol': official_gene_symbol, 'official protein name': official_protein_name, 
         'uniprot AC': uniprot_AC, 'GeneID': gene_ids, 'description': description}
    Basic_info = pd.DataFrame(data=d)        
    Basic_info.to_csv('Basic_Information.csv', sep='\t')
    
    return

