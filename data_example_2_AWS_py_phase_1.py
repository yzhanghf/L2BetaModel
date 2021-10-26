import sys
import os
import time
import csv

import numpy as np
from numpy import linalg as LA

import pandas as pd

node_list_df = pd.read_csv('data/AWS-csv/nodes/paper_nodes.csv')
edge_list_df = pd.read_csv('data/AWS-csv/edges/paper_to_reference.csv')

node_name_list = node_list_df['~id']
paper_dict = dict.fromkeys(node_name_list,0)
for iii in range(len(paper_dict)):
    paper_dict[node_name_list[iii]] = iii

from_list = edge_list_df['~from']
to_list   = edge_list_df['~to']

deg_dict = dict.fromkeys(node_name_list,0)
for iii in range(len(from_list)):
    deg_dict[from_list[iii]] += 1
    deg_dict[to_list[iii]]   += 1

deg_vec = list(deg_dict.values())

with open('data/AWS-csv/deg_vec.txt','w') as wt:
    for ii in range(N):
        _ = wt.write('%d\n'%(deg_vec[ii]));

with open('data/AWS-csv/node_name_vec.txt','w') as wt:
    for ii in range(N):
        _ = wt.write('%s\n'%(node_name_list[ii]));

#################################
# construct institution -> country dictionary

inst_to_natl_df = pd.read_csv('data/AWS-csv/nodes/institution_nodes.csv')

nation_list  = inst_to_natl_df['country:String']
id_list = inst_to_natl_df['~id']
# NameList = ['China','US','UK','EU','JpKr','India','beta_est']

institution_dict = dict.fromkeys(id_list, 0)

for iii in range(len(institution_dict)):
    temp = nation_list[iii]
    
    if type(temp) is not str:
        continue
    
    if 'China' in temp or 'china' in temp:
        institution_dict[id_list[iii]] = 1
    
    if 'U' in temp and 'S' in temp:
        institution_dict[id_list[iii]] = 2
    
    if 'U' in temp and 'K' in temp:
        institution_dict[id_list[iii]] = 3
    
    if 'France' in temp or 'Germany' in temp or 'Italy' in temp or 'Spain' in temp or 'Netherlands' in temp or 'Switzerland' in temp:
        institution_dict[id_list[iii]] = 4
    
    if 'Japan' in temp or 'Korea' in temp:
        institution_dict[id_list[iii]] = 5
    
    if 'India' in temp:
        institution_dict[id_list[iii]] = 6

# construct author -> country dictionary
author_list_df = pd.read_csv('data/AWS-csv/nodes/paper_author_nodes.csv')
author_list = author_list_df['~id']
author_list_dict = dict.fromkeys(author_list,0)
for iii in range(len(author_list)):
    author_list_dict[author_list[iii]]=iii

author_nation_table = np.zeros([len(author_list),7])


author_to_institution_df = pd.read_csv('data/AWS-csv/edges/author_to_institution.csv')
x1 = author_to_institution_df['~from'] # authors
x2 = author_to_institution_df['~to'] # institutions

for iii in range(len(author_to_institution_df)):
    author_nation_table[author_list_dict[x1[iii]],institution_dict[x2[iii]]] = 1


# construct paper -> nation table
paper_nation_table = np.zeros([len(paper_dict), 7])


paper_to_author_df = pd.read_csv('data/AWS-csv/edges/paper_to_author.csv')
y1 = paper_to_author_df['~from'] # papers
y2 = paper_to_author_df['~to'] # authors

for iii in range(len(y1)):
    paper_nation_table[paper_dict[y1[iii]], :] += author_nation_table[author_list_dict[y2[iii]], :]


paper_nation_table = 1*(paper_nation_table>0);

np.savetxt('./data/AWS-csv/nation.csv', paper_nation_table, delimiter=',')
























