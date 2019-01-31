#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import Interactions
import pandas as pd
import networkx as nx
import markov_clustering as mc
import community
from scipy.stats import hypergeom
import os



def adjacency_matrix(SG, biogrid, iid):
    print('Netwrok Analysing...')
    sgi, union_ = Interactions._integrations(SG, biogrid, iid)
    intersections = Interactions.intersection_interactions(SG, biogrid, iid)
    
    # making adjacency matrix of union_integrated dataset
    union_adjacency = pd.crosstab(union_['interactor A'], union_['interactor B'])
    idx = union_adjacency.columns.union(union_adjacency.index)
    union_adjacency = union_adjacency.reindex(index = idx, columns=idx, fill_value=0)
    
    # making adjacency matrix of sgi_integrated dataset
    sgi_adjacency = pd.crosstab(sgi['interactor A'], sgi['interactor B'])
    idx = sgi_adjacency.columns.union(sgi_adjacency.index)
    sgi_adjacency = sgi_adjacency.reindex(index = idx, columns=idx, fill_value=0)
    
    # making adjacency matrix of intersections dataset
    I_adjacency = pd.crosstab(intersections['interactor A'], intersections['interactor B'])
    idx = I_adjacency.columns.union(I_adjacency.index)
    I_adjacency = I_adjacency.reindex(index = idx, columns=idx, fill_value=0)
    print('adjacency matrices has been made!')
    return union_adjacency, sgi_adjacency, I_adjacency

def create_graph(union_adjacency, sgi_adjacency, I_adjacency):
    
    # making complete graph of union data_set
    G_union = nx.from_pandas_adjacency(union_adjacency)
    G_union.name = 'Union_Interactions Graph'
    
    # making complete graph of sgi data_set
    G_sgi = nx.from_pandas_adjacency(sgi_adjacency)
    G_sgi.name = 'SGI Graph'
    
    # making complete graph of Intersection data_set
    G_I = nx.from_pandas_adjacency(I_adjacency)
    G_I.name = 'Intersection interactions Graph'
    
    # making large connected component of union graph
    G_lcc_union = max(nx.connected_component_subgraphs(G_union), key=len)
    G_lcc_union.name = 'Large Connected Component of G_union'
    
    G_lcc_I = max(nx.connected_component_subgraphs(G_I), key=len)
    G_lcc_I.name = 'Large Connected Component of G_I'
    print('Graphs Created!')
    return G_union, G_sgi, G_I, G_lcc_union, G_lcc_I


def mcl_eval(cluster, list_of_interactomes, SG):
    # Finding seed_genes blongs to each cluster of MCL of union set
    # cluster is a list of tuples representing each cluster wich contains index of each node in the main graph
    seed_genes_belongs_to_each_cluster = [[] for x in range(len(cluster))]
    for i in range(len(cluster)):
        for j in range(len(cluster[i])):
            if(list_of_interactomes[cluster[i][j]] in SG):
                # by taking the index of related Uniport in SG_uniprot list, since the length of uniprot list and SG list is the same
                # we just append Seed_Gene symbol in the final list, cause seed_gene symbol is more important and understandable
                seed_genes_belongs_to_each_cluster[i].append(list_of_interactomes[cluster[i][j]])
                
    number_of_seed_genes_in_each_cluster = [len(l) for l in seed_genes_belongs_to_each_cluster]
    
    # total number of interactomes in each cluster        
    clusters_length = [len(t) for t in cluster]
        
    return seed_genes_belongs_to_each_cluster, number_of_seed_genes_in_each_cluster, clusters_length

def louvain_eval(partition, SG):
    
    # communities list is a list of lists that represents the uniprots blongs to each partition
    communities_list = [[] for x in range(max(partition.values()))]
    for i in range(max(partition.values())):
        for n,p in partition.items():
            if (p == i):
                communities_list[i].append(n)
    SG_of_partitions = [[] for x in range(len(communities_list))]
    for i in range(len(communities_list)):
        for j in range(len(communities_list[i])):
            if (communities_list[i][j] in SG):
                SG_of_partitions[i].append(communities_list[i][j])
    
    number_of_seed_genes_in_each_partition = [len(l) for l in SG_of_partitions]
    
    partitions_length = [len(l) for l in communities_list]
    
    return len(communities_list), SG_of_partitions, number_of_seed_genes_in_each_partition, partitions_length
    
    
def clustering(G_lcc_union, G_lcc_I, SG):
    # MCL clustering Union_lcc 
    matrix_union = nx.to_scipy_sparse_matrix(G_lcc_union)
    result_union_mcl = mc.run_mcl(matrix_union)           
    clusters_union = mc.get_clusters(result_union_mcl)
    
    # MCL clustering Intersection_lcc 
    matrix_I = nx.to_scipy_sparse_matrix(G_lcc_I)
    result_I_mcl = mc.run_mcl(matrix_I)           
    clusters_I = mc.get_clusters(result_I_mcl)
    
    # Louvain clustering Union_lcc
    partition_u = community.best_partition(G_lcc_union)
    # Louvain clustering Intersection_lcc
    partition_I = community.best_partition(G_lcc_I)
    
    
    # MCL_eval
    list_of_genes_u = list(G_lcc_union.nodes())
    SG_of_each_cluster_u, num_SG_of_cluster_u, clusters_length_u = mcl_eval(clusters_union, list_of_genes_u, SG)
    list_of_genes_I = list(G_lcc_I.nodes())
    SG_of_each_cluster_I, num_SG_of_cluster_I, clusters_length_I = mcl_eval(clusters_I, list_of_genes_I, SG)
    # Louvain_eval
    num_of_communities_u, SG_of_partitions_u, num_SG_of_partitions_u, partitions_length_u = louvain_eval(partition_u, SG)
    num_of_communities_I, SG_of_partitions_I, num_SG_of_partitions_I, partitions_length_I = louvain_eval(partition_I, SG)
    
    
    M_u = 80822 #total no. of union interactions
    M_I = 914 #total no. of intersection interactions
    n = 27 # no. of seed_genes
    
    # preparing results in pandas format
    module_number_mcl_u = [i for i in range(len(clusters_union))]
    pval_mcl_u = [hypergeom.sf(num_SG_of_cluster_u[i]-1, M_u, n, len(clusters_union[i])) for i in range(len(clusters_union))]
    ratio_mcl_u = [round(num_SG_of_cluster_u[i]/clusters_length_u[i], 4) for i in range(len(clusters_union))]
    putative_mcl_u = ['Yes' if p < 0.05 else 'No' for p in pval_mcl_u]
    d_mcl_u = {'Module Number': module_number_mcl_u, 'Number of Genes': clusters_length_u, 'Number of Seed_genes': num_SG_of_cluster_u,
               'Ratio': ratio_mcl_u, 'p-Value': pval_mcl_u, 'Putative disease': putative_mcl_u}
    pd_mcl_u = pd.DataFrame(data = d_mcl_u)
    
    module_number_mcl_I = [i for i in range(len(clusters_I))]
    # calculating p-val of each cluster using hypergeom library
    pval_mcl_I = [hypergeom.sf(num_SG_of_cluster_I[i]-1, M_I, n, len(clusters_I[i])) for i in range(len(clusters_I))]
    # ration (in each cluster) = number of seed genes / total number of genes 
    ratio_mcl_I = [round(num_SG_of_cluster_I[i]/clusters_length_I[i], 4) for i in range(len(clusters_I))]
    # putative disease refers to the cluster with significant p-val (pval<0.05)
    putative_mcl_I = ['Yes' if p < 0.05 else 'No' for p in pval_mcl_I]
    d_mcl_I = {'Module Number': module_number_mcl_I, 'Number of Genes': clusters_length_I, 'Number of Seed_genes': num_SG_of_cluster_I,
               'Ratio': ratio_mcl_I, 'p-Value': pval_mcl_I, 'Putative disease': putative_mcl_I}
    pd_mcl_I = pd.DataFrame(data = d_mcl_I)
    
    module_number_louvain_u = [i for i in range(num_of_communities_u)]
    pval_louvain_u = [hypergeom.sf(num_SG_of_partitions_u[i]-1, M_u, n, partitions_length_u[i]) for i in range(num_of_communities_u)]
    ratio_louvain_u = [round(num_SG_of_partitions_u[i]/partitions_length_u[i], 4) for i in range(num_of_communities_u)]
    putative_louvain_u = ['Yes' if p < 0.05 else 'No' for p in pval_louvain_u]
    d_louvain_u = {'Module Number': module_number_louvain_u, 'Number of Genes': partitions_length_u, 'Number of Seed_genes': num_SG_of_partitions_u,
               'Ratio': ratio_louvain_u, 'p-Value': pval_louvain_u, 'Putative disease': putative_louvain_u}
    pd_louvain_u = pd.DataFrame(data = d_louvain_u)
    
    module_number_louvain_I = [i for i in range(num_of_communities_I)]
    pval_louvain_I = [hypergeom.sf(num_SG_of_partitions_I[i]-1, M_I, n, partitions_length_I[i]) for i in range(num_of_communities_I)]
    ratio_louvain_I = [round(num_SG_of_partitions_I[i]/partitions_length_I[i], 4) for i in range(num_of_communities_I)]
    putative_louvain_I = ['Yes' if p < 0.05 else 'No' for p in pval_louvain_I]
    d_louvain_I = {'Module Number': module_number_louvain_I, 'Number of Genes': partitions_length_I, 'Number of Seed_genes': num_SG_of_partitions_I,
               'Ratio': ratio_louvain_I, 'p-Value': pval_louvain_I, 'Putative disease': putative_louvain_I}
    pd_louvain_I = pd.DataFrame(data = d_louvain_I)
    
    print('Clustering has been done!')
    return pd_mcl_u, pd_mcl_I, pd_louvain_u, pd_louvain_I
    
def save_results(G_union, G_sgi, G_I, G_lcc_union, G_lcc_I, pd_mcl_u, pd_mcl_I, pd_louvain_u, pd_louvain_I):

    os.mkdir('Network_Analysis_results', mode = 0o777)
    os.chdir('Network_Analysis_results')
    
    f = open('Network_meesures.txt', 'w')
    
    ### Results of G_Union
    f.write('Union_complete_graph details....................\n')
    f.write(nx.info(G_union))
    f.write('No. of connected component: '+ str(nx.number_connected_components(G_union))+ '\n')
    f.write('No. of isolated Nodes: '+ str(len(list(nx.isolates(G_union))))+ '\n')
    avg_shortest_path_u = [nx.average_shortest_path_length(g) for g in nx.connected_component_subgraphs(G_union)]
    f.write('Average_path_length: '+ '\n')
    for i in range(len(avg_shortest_path_u)):
        f.write(str(avg_shortest_path_u[i])+ '\n')
    
    # The average degree connectivity is the average nearest neighbor degree of nodes with degree k. 
    avg_deg_union = nx.average_degree_connectivity(G_union)
    avg_deg_union_res = open('average_degree_union.txt', 'w')
    avg_deg_union_res.write('\t'.join(['#degree k','average connectivity'])+ '\n')
    for k, v in avg_deg_union.items():
        avg_deg_union_res.write('\t'.join(map(str,([k,v])))+ '\n')
    avg_deg_union_res.close()
    
    f.write('Average clustering coefficient: '+ str(nx.average_clustering(G_union))+ '\n')
    # The diameter is the maximum eccentricity.
    diameter_u = [nx.diameter(g) for g in nx.connected_component_subgraphs(G_union)]
    f.write('Network diameter: ' + '\n')
    for i in range(len(diameter_u)):
        f.write(str(diameter_u[i])+ '\n')
    # The radius is the minimum eccentricity.
    radius_u = [nx.radius(g) for g in nx.connected_component_subgraphs(G_union)]
    f.write('Network radius: ' + '\n')
    for i in range(len(radius_u)):
        f.write(str(radius_u[i]) + '\n')
    
    
    ### Results of G_sgi
    f.write('SGI_complete_graph details....................\n')
    f.write(nx.info(G_sgi))
    f.write('No. of connected component: '+ str(nx.number_connected_components(G_sgi))+ '\n')
    f.write('No. of isolated Nodes: '+ str(len(list(nx.isolates(G_sgi))))+ '\n')
    avg_shortest_path_sgi = [nx.average_shortest_path_length(g) for g in nx.connected_component_subgraphs(G_sgi)]
    f.write('Average_path_length: '+ '\n')
    for i in range(len(avg_shortest_path_sgi)):
        f.write(str(avg_shortest_path_sgi[i]) + '\n')
    
    # The average degree connectivity is the average nearest neighbor degree of nodes with degree k. 
    avg_deg_sgi = nx.average_degree_connectivity(G_sgi)
    avg_deg_sgi_res = open('average_degree_sgi.txt', 'w')
    avg_deg_sgi_res.write('\t'.join(['#degree k','average connectivity'])+ '\n')
    for k, v in avg_deg_sgi.items():
        avg_deg_sgi_res.write('\t'.join(map(str,([k,v])))+ '\n')
    avg_deg_sgi_res.close()
    
    f.write('Average clustering coefficient: '+ str(nx.average_clustering(G_sgi))+ '\n')
    # The diameter is the maximum eccentricity.
    diameter_sgi = [nx.diameter(g) for g in nx.connected_component_subgraphs(G_sgi)]
    f.write('Network diameter: ' + '\n')
    for i in range(len(diameter_sgi)):
        f.write(str(diameter_sgi[i])+ '\n')
    # The radius is the minimum eccentricity.
    radius_sgi = [nx.radius(g) for g in nx.connected_component_subgraphs(G_sgi)]
    f.write('Network radius: ' + '\n')
    for i in range(len(radius_sgi)):
        f.write(str(radius_sgi[i]) + '\n')
    
    ### Results of G_I
    f.write('Intersections_complete_graph details....................\n')
    f.write(nx.info(G_I))
    f.write('No. of connected component: '+ str(nx.number_connected_components(G_I))+ '\n')
    f.write('No. of isolated Nodes: '+ str(len(list(nx.isolates(G_I))))+ '\n')
    avg_shortest_path_I = [nx.average_shortest_path_length(g) for g in nx.connected_component_subgraphs(G_I)]
    f.write('Average_path_lengths: '+ '\n')
    for i in range(len(avg_shortest_path_I)):
        f.write(str(avg_shortest_path_I[i])+ '\n')
    
    # The average degree connectivity is the average nearest neighbor degree of nodes with degree k. 
    avg_deg_I = nx.average_degree_connectivity(G_I)
    avg_deg_I_res = open('average_degree_Intersection.txt', 'w')
    avg_deg_I_res.write('\t'.join(['#degree k','average connectivity'])+ '\n')
    for k, v in avg_deg_I.items():
        avg_deg_I_res.write('\t'.join(map(str,([k,v])))+ '\n')
    avg_deg_I_res.close()
    
    f.write('Average clustering coefficient: '+ str(nx.average_clustering(G_I))+ '\n')
    # The diameter is the maximum eccentricity.
    diameter_I = [nx.diameter(g) for g in nx.connected_component_subgraphs(G_I)]
    f.write('Network diameter: ' + '\n')
    for i in range(len(diameter_I)):
        f.write(str(diameter_I[i])+ '\n')
    # The radius is the minimum eccentricity.
    radius_I = [nx.radius(g) for g in nx.connected_component_subgraphs(G_I)]
    f.write('Network radius: ' + '\n')
    for i in range(len(radius_I)):
        f.write(str(radius_I[i]) + '\n')
    
    f.close()
    
    
    f_lcc = open('lcc_Network_measures.txt', 'w')
    f_lcc.write('#### Large Connected Component Graphs ####\n')
    
    # Results of G_lcc_union
    f_lcc.write('Union_lcc_complete_graph details....\n')
    f_lcc.write(nx.info(G_lcc_union))
    f_lcc.write('Average_path_length: '+ str(nx.average_shortest_path_length(G_lcc_union))+ '\n')
    
    # The average degree connectivity is the average nearest neighbor degree of nodes with degree k. 
    avg_deg_union_lcc = nx.average_degree_connectivity(G_lcc_union)
    avg_deg_unionlcc_res = open('average_degree_Union_lcc.txt', 'w')
    avg_deg_unionlcc_res.write('\t'.join(['#degree k','average connectivity'])+ '\n')
    for k, v in avg_deg_union_lcc.items():
        avg_deg_unionlcc_res.write('\t'.join(map(str,([k,v])))+ '\n')
    avg_deg_unionlcc_res.close()
    
    f_lcc.write('Average clustering coefficient: '+ str(nx.average_clustering(G_lcc_union))+ '\n')
    # The diameter is the maximum eccentricity.
    f_lcc.write('Network diameter: ' + str(nx.diameter(G_lcc_union))+ '\n')
    # The radius is the minimum eccentricity.
    f_lcc.write('Network radius: ' + str(nx.radius(G_lcc_union)) + '\n')
    
    # Results of G_lcc_I
    f_lcc.write('Intersection_lcc_complete_graph details....\n')
    f_lcc.write(nx.info(G_lcc_I))
    f_lcc.write('Average_path_length: '+ str(nx.average_shortest_path_length(G_lcc_I))+ '\n')
    
    # The average degree connectivity is the average nearest neighbor degree of nodes with degree k. 
    avg_deg_I_lcc = nx.average_degree_connectivity(G_lcc_I)
    avg_deg_I_lcc_res = open('average_degree_Intersection_lcc.txt', 'w')
    avg_deg_I_lcc_res.write('\t'.join(['#degree k','average connectivity'])+ '\n')
    for k, v in avg_deg_I_lcc.items():
        avg_deg_I_lcc_res.write('\t'.join(map(str,([k,v])))+ '\n')
    avg_deg_I_lcc_res.close()
    
    f_lcc.write('Average clustering coefficient: '+ str(nx.average_clustering(G_lcc_I))+ '\n')
    # The diameter is the maximum eccentricity.
    f_lcc.write('Network diameter: ' + str(nx.diameter(G_lcc_I))+ '\n')
    # The radius is the minimum eccentricity.
    f_lcc.write('Network radius: ' + str(nx.radius(G_lcc_I)) + '\n')
    
    f.close()
    
    pd_mcl_u.to_csv('MCL_union_results.csv', sep='\t')
    pd_mcl_I.to_csv('MCL_intersection_results.csv', sep='\t')
    pd_louvain_u.to_csv('louvain_union_results.csv', sep='\t') 
    pd_louvain_I.to_csv('louvain_intersection_results.csv', sep='\t')
    os.chdir('..')
    print('Network Analysing has been done successfully and results saved in a proper directory!')
    return
    