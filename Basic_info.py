#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import pandas as pd
import urllib.request as ur
from bs4 import BeautifulSoup as bs

# obtaining Basic informations from web and local datasets
def Basic_informations(SG, iid):
    gene_ids = list()
    uniprot_AC = list()
    official_gene_symbol = list()
    official_protein_name = list()
    description = list()
    print('Fetching Basic Informations...')
    # the first loop is to get the gene IDs from ncbi.nlm.nih.gov which is also approved by HGNC
    for i in range(len(SG)):
        ids = list()
        em_list = list()
        index_list = list()
        response = ur.urlopen('https://www.ncbi.nlm.nih.gov/gene/?term=' + SG[i])
        html = response.read()
        soup = bs(html, 'html.parser')
        for items in soup.findAll('em'):
            em_list.append(items.text[0:4])
            for i in range(len(em_list)):
                if (em_list[i] == 'Homo'):
                    index_list.append(i)
        for items in soup.findAll('span', {'class', 'gene-id'}):
            ids.append(items.text[4:])
        gene_ids.append(ids[index_list[0]])
    # second loop is to get the proper Uniprot_AC of each seed gene by processing local iid data set 
    for i in range(len(SG)):
        uniprot_AC.append(iid[iid.symbol1 == SG[i]].iloc[0][2])
    # third (and final) loop is to go through all the fetched gene_ids in first part (open the proper webpage on ncbi.nlm.nih.gov) and fetch the interested
    # informations: official_gene_symbol, official_protein_name, breif_description
    for i in range(len(gene_ids)):
        desc_list = list()
        response = ur.urlopen('https://www.ncbi.nlm.nih.gov/gene/'+ gene_ids[i])
        html = response.read()
        soup = bs(html, 'html.parser')
        for items in soup.findAll(['dd']):
            desc_list.append(items.text)
        official_gene_symbol.append(desc_list[0][0:-16])
        official_protein_name.append(desc_list[1][0:-16])
        description.append(desc_list[9].split('.')[0])                  
    print('Basic Informations has been fetched successfully!')   
    d = {'seed_genes': SG,'official gene symbol': official_gene_symbol, 'official protein name': official_protein_name, 
         'uniprot AC': uniprot_AC, 'GeneID': gene_ids, 'description': description}
    Basic_info = pd.DataFrame(data=d)        
    Basic_info.to_csv('Basic_Information.csv', sep='\t')
    
    