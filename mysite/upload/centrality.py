

__author__ = "Zhaolong Yu"
__copyright__ = "Copyright 2017"
__credits__ = ["Zhaolong Yu"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Zhaolong Yu"
__email__ = "zhaolong.yu@yale.edu"

### Usage:      python3 centrality_calculation.py -i <input file> -a <annotation file> -n <dip ppi file>
### Example:    python3 centrality_calculation.py -i input.txt -a map_table.txt -n dip.txt


import argparse
import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
import networkx as nx


def calculate(args_input):
    ## open the PPI file
    fo_dip = open('media/datasets/dip.txt',"r")
    fo_dip.readline()
    edges = []
    dip_uniprot = {}
    uniprot_dip = {}
    for rawline in fo_dip:
        line = rawline.strip().split("\t")
        list1 = line[0].split("|")
        list2 = line[1].split("|")
        if len(list1)>0 and len(list2)>0:
            node1 = line[0].split("|")[0]
            node2 = line[1].split("|")[0]
            edges.append([node1,node2])
        ref = list1[-1].split(":")
        if len(list1)>1 and ref[0] == "uniprotkb":
            if list1[0] not in dip_uniprot.keys():
                dip_uniprot[list1[0]] = ref[1]
                uniprot_dip[ref[1]] = list1[0]
        ref = list2[-1].split(":")
        if len(list2)>1 and ref[0] == "uniprotkb":
            if list2[0] not in dip_uniprot.keys():
                dip_uniprot[list2[0]] = ref[1]
                uniprot_dip[ref[1]] = list2[0]
    
            
        
    fo_dip.close()
    
    ## genes in the PPI network
    all_genes = []
    for gene in dip_uniprot.items():
        all_genes.append(gene[1])
    
        
    ## read the snp file
    #fo_snp = open("Z.3DStruct_annotation.txt","r")
    fo_snp = open(args_input,"r")
    
    snp_genes = []
    for rawline in fo_snp:
        line = rawline.strip().split("\t")
        gene = line[2].split(";")[1]
        snp_genes.append(gene)
    fo_snp.close()   
    
    ## get the proteins containing SNPs 
    unique_snp_genes = np.unique(snp_genes)
      
    ## gene annotation file
    #fo_map = open("map_table.txt","r")
    fo_map = open("media/datasets/map_table.txt","r")
    fo_map.readline()
    gene_uniprot = {}
    uniprot_gene = {}
    for rawline in fo_map:
        line = rawline.strip().split("\t")
        gene_uniprot[line[1]] = line[0]
        uniprot_gene[line[0]] = line[1]
        
    fo_map.close()
    
    ## get the gene list
    gene_list = []
    for gene in all_genes:
        if gene in uniprot_gene.keys():
            gene_list.append(uniprot_gene[gene])
    
    ## gene set partitioning
    unique_nonsnp_genes = set(gene_list)-set(unique_snp_genes)       
    unique_snp_genes = set(unique_snp_genes).intersection(set(gene_list))        
        
    
    ## building the network
    g = nx.Graph() 
    g.add_edges_from(edges)   
    
    g.number_of_nodes()  
    
    snp_degree_centrality = {}
    nonsnp_degree_centrality = {}
    
    snp_betweenness_centrality = {}
    nonsnp_betweenness_centrality = {}
    
    degree_dict = {}
    degree_dict = nx.degree(g)
        
    ## calculate the degree centrality
    tmp_dc = nx.degree_centrality(g)
    
    ## result dict for degree centrality
    for gene in unique_nonsnp_genes:
        nonsnp_degree_centrality[gene] = tmp_dc[uniprot_dip[gene_uniprot[gene]]]
    
    for gene in unique_snp_genes:
        snp_degree_centrality[gene] = tmp_dc[uniprot_dip[gene_uniprot[gene]]]
    
    ## calculate the betweenness centrality
    tmp_bc = nx.betweenness_centrality(g,k=4901,normalized=True)
    
    ## result dict for betweenness centrality    
    for gene in unique_nonsnp_genes:
        nonsnp_betweenness_centrality[gene] = tmp_bc[uniprot_dip[gene_uniprot[gene]]]
         
    for gene in unique_snp_genes:
        snp_betweenness_centrality[gene] = tmp_bc[uniprot_dip[gene_uniprot[gene]]]
    
    ## result list (for further analysis)    
    snp_degree_centrality_list = []
    nonsnp_degree_centrality_list = []
    
    snp_betweenness_centrality_list = []
    nonsnp_betweenness_centrality_list = []
          
    for gene in snp_degree_centrality.keys():
        snp_degree_centrality_list.append(float(snp_degree_centrality[gene]))
        
    for gene in nonsnp_degree_centrality.keys():
        nonsnp_degree_centrality_list.append(float(nonsnp_degree_centrality[gene]))
        
    for gene in snp_betweenness_centrality.keys():
        snp_betweenness_centrality_list.append(float(snp_betweenness_centrality[gene]))
        
    for gene in nonsnp_betweenness_centrality.keys():
        nonsnp_betweenness_centrality_list.append(float(nonsnp_betweenness_centrality[gene]))
    
    ## sorting results and output the top10 list
    #f = open('media/datasets/nonsnp_centrality.csv','w')
    f = open('personal/static/datasets/nonsnp_centrality.csv','w')
    f.write("Gene" + "," + "Uniprot_ID" + "," + "DIP_ID" + "," + "Degree_Centrality" + "," + "Betweenness_Centrality" + "\n")
    sorted_keys = sorted(nonsnp_degree_centrality, key=nonsnp_degree_centrality.get, reverse=True)
    for r in sorted_keys:
        f.write(r + "," + gene_uniprot[r] + "," + uniprot_dip[gene_uniprot[r]] + "," + str(nonsnp_degree_centrality[r])+ ","  +str(nonsnp_betweenness_centrality[r]) +  "\n")
    f.close()

    f = open('personal/static/datasets/snp_centrality.csv','w')
    f.write("Gene" + "," + "Uniprot_ID" + "," + "DIP_ID" + "," + "Degree_Centrality" + "," + "Betweenness_Centrality" + "\n")
    sorted_keys = sorted(snp_degree_centrality, key=snp_degree_centrality.get, reverse=True)
    for r in sorted_keys:
        f.write(r + "," + gene_uniprot[r] + "," + uniprot_dip[gene_uniprot[r]] + "," + str(snp_degree_centrality[r]) + "," + str(snp_betweenness_centrality[r]) +  "\n")
    f.close()

