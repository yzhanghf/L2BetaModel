import sys
import os
import time
import csv

import numpy as np
from numpy import linalg as LA

import pandas as pd

import networkx as nx

M = 10179046

kg_dict = dict()
kg_act  = dict()
act_vec = np.zeros([M,1])
line_list = np.zeros([M,2])

line = 123;count = 0
fp = open('data/kg.nt.original', 'r', encoding='utf-8')
while line:
    
    if count % 100000 == 0:
        print(count, end=' ', flush=True);
    
    line = fp.readline()
    x = line.split('> <')
    
    if len(x)<=2:
        continue
    
    temp1 = x[0]+'>'
    temp_act = x[1];
    temp2 = '<'+x[2]
    temp2 = temp2[:-3]
    
    if temp1[0] != '<' or temp2[-1] != '>':
        break
    
    if temp1 == temp2:
        continue
    
    temp_act = x[1];
    if temp_act not in kg_act:
        x0 = len(kg_act)
        kg_act[temp_act] = x0
        act_vec[count] = x0
    else:
        act_vec[count] = kg_act[temp_act]
    
    if temp1 not in kg_dict:
        x1 = len(kg_dict)
        kg_dict[temp1] = x1
    else:
        x1 = kg_dict[temp1]
        
    if temp2 not in kg_dict:
        x2 = len(kg_dict)
        kg_dict[temp2] = x2
    else:
        x2 = kg_dict[temp2]
    
    line_list[count,:] = [x1,x2]
    
    count += 1

fp.close()

act_vec = act_vec[range(count)];
line_list = line_list[range(count),:];

################################################
################################################

keylist = list(kg_dict.keys())


################################################
################################################

zz00 = (np.round(line_list[np.where(act_vec==0)[0],0])).astype(int)
zz01 = (np.round(line_list[np.where(act_vec==0)[0],1])).astype(int)

zz2 = [keylist[iii] for iii in list(set(zz01))]

_ = [print(kkk) for kkk in zz2]

_ = [print(kkk) for kkk in kg_act.keys()]

for kkk in range(len(zz2)):
    temp_idx = [jjj for jjj in range(len(zz01)) if keylist[zz01[jjj]]==zz2[kkk]]
    for iii in range(min([5,len(temp_idx)])):
        print(keylist[zz00[temp_idx[iii]]], end='  ')
        print(keylist[zz01[temp_idx[iii]]])
    print('\n')

zz3 = ['Work', 'JournalArticle', 'Journal', 'Book', 'ReviewPaper', 'BookSeries', 'BookChapter', 'Organization', 'Person', 'Thesis', 'ReferenceWork', 'BookSet', 'ProceedingsPaper', 'Dataset', 'ConferenceProceedings']

zz_name_dict = dict();
for iii in range(len(zz2)):
    zz_name_dict[zz2[iii]] = zz3[iii]



################################################
################################################

# creat paper list

paper_dict  = dict()
person_dict = dict()
org_dict    = dict()


line = 123;count = 0
fp = open('data/kg.nt.original', 'r', encoding='utf-8')
while line:
    
    if count % 100000 == 0:
        print(count, end=' ', flush=True);
    
    line = fp.readline()
    x = line.split('> <')
    
    if len(x)<=2:
        continue
    
    temp1 = x[0]+'>'
    temp_act = x[1];
    temp2 = '<'+x[2]
    temp2 = temp2[:-3]
    
    if temp1[0] != '<' or temp2[-1] != '>':
        break
    
    if temp1 == temp2:
        continue
    
    if 'dx.doi.org' in temp1:
        if temp1 not in paper_dict:
            paper_dict[temp1] = len(paper_dict)
    
    if 'dx.doi.org' in temp2:
        if temp2 not in paper_dict:
            paper_dict[temp2] = len(paper_dict)
    
    
    if 'orcid.org' in temp1:
        if temp1 not in person_dict:
            person_dict[temp1] = len(person_dict)
    
    if 'orcid.org' in temp2:
        if temp2 not in person_dict:
            person_dict[temp2] = len(person_dict)
    
    
    if 'idlab.github.io' in temp1:
        if temp1 not in org_dict:
            org_dict[temp1] = len(org_dict)
    
    if 'idlab.github.io' in temp2:
        if temp2 not in org_dict:
            org_dict[temp2] = len(org_dict)
    
    count += 1

fp.close()




# Construct degree list
N = len(paper_dict)
deg_vec = np.zeros([N,1], dtype=int)

line = 123;count = 0
fp = open('data/kg.nt.original', 'r', encoding='utf-8')
while line:
    
    if count % 100000 == 0:
        print(count, end=' ', flush=True);
    
    line = fp.readline()
    x = line.split('> <')
    
    if len(x)<=2:
        continue
    
    temp1 = x[0]+'>'
    temp_act = x[1];
    temp2 = '<'+x[2]
    temp2 = temp2[:-3]
    
    if temp1[0] != '<' or temp2[-1] != '>':
        break
    
    if temp1 == temp2:
        continue
    
    if 'dx.doi.org' in temp1 and 'dx.doi.org' in temp2:
        if 'cite' not in temp_act:
            print('Error!\n')
        
        deg_vec[paper_dict[temp1]] += 1
        deg_vec[paper_dict[temp2]] += 1
    
    count += 1

fp.close()

with open('data/kg_degs.txt','w') as wt:
    for ii in range(N):
        _ = wt.write('%d\n'%(deg_vec[ii]))

paper_keylist = list(paper_dict.keys())
with open('data/kg_dict_keys.txt','w') as wt:
    for ii in range(N):
        _ = wt.write('%s\n'%(paper_keylist[ii]))


## NEXT TO DO: CONSTRUCT NODAL COVARIATES

person_to_org_dict  = dict()
org_to_city_dict    = dict()
city_to_nation_dict = dict()
org_list    = list()
city_list   = list()
nation_list = list()

line = 123;count = 0
fp = open('data/kg.nt.original', 'r', encoding='utf-8')
while line:
    
    if count % 100000 == 0:
        print(count, end=' ', flush=True);
    
    line = fp.readline()
    x = line.split('> <')
    
    if len(x)<=2:
        continue
    
    temp1 = x[0]+'>'
    temp_act = x[1];
    temp2 = '<'+x[2]
    temp2 = temp2[:-3]
    
    if kg_act[temp_act] == 5 and temp2 != '<http://idlab.github.io/covid19#>': # author -> institute(org)
        if temp1 not in person_to_org_dict:
            person_to_org_dict[temp1] = temp2
        else:
            if type(person_to_org_dict[temp1]) is not list:
                if temp2 != person_to_org_dict[temp1]:
                    person_to_org_dict[temp1] = list([person_to_org_dict[temp1], temp2])
            else:
                if temp2 not in person_to_org_dict[temp1]:
                    person_to_org_dict[temp1].append(temp2)
        
        if temp2 not in org_list:
            org_list.append(temp2)
        
        
    if kg_act[temp_act] == 3:
        if temp1 not in org_to_city_dict:
            org_to_city_dict[temp1] = temp2
        else:
            if type(org_to_city_dict[temp1]) is not list:
                if temp2 != org_to_city_dict[temp1]:
                    org_to_city_dict[temp1] = list([org_to_city_dict[temp1], temp2])
            else:
                if temp2 not in org_to_city_dict[temp1]:
                    org_to_city_dict[temp1].append(temp2)
        
        if temp1 not in org_list:
            org_list.append(temp1)
        
        if temp2 not in city_list:
            city_list.append(temp2)
    
    if kg_act[temp_act] == 4:
        if temp1 not in city_to_nation_dict:
            city_to_nation_dict[temp1] = temp2
        else:
            if type(city_to_nation_dict[temp1]) is not list:
                if temp2 != city_to_nation_dict[temp1]:
                    city_to_nation_dict[temp1] = list([city_to_nation_dict[temp1], temp2])
            else:
                if temp2 not in city_to_nation_dict[temp1]:
                    city_to_nation_dict[temp1].append(temp2)
        
        if temp1 not in city_list:
            city_list.append(temp1)
        
        if temp2 not in nation_list:
            nation_list.append(temp2)
    
    
    count += 1

fp.close()



############################################


All_Nodes = pd.DataFrame(columns = ['DOI','degree','type','China','USA','UK','EU','JapanKorea','India'])
All_Nodes.DOI = paper_dict.keys()
All_Nodes.degree = deg_vec
All_Nodes['type'] = '';



def Func_nation_count(string_var):
    count_vec = np.zeros([1,6], dtype=int)
    if type(string_var) is not list:
        if 'China' in string_var:
            count_vec[0,0] += 1
        
        if 'U' in string_var and 'S' in string_var and 'A' in string_var:
            count_vec[0,1] += 1
        
        if 'U' in string_var and 'K' in string_var:
            count_vec[0,2] += 1
        
        if 'Germany' in string_var or 'France' in string_var or 'Italy' in string_var or 'Spain' in string_var:
            count_vec[0,3] += 1
        
        if 'Japan' in string_var or 'Korea' in string_var:
            count_vec[0,4] += 1
        
        if 'India' in string_var:
            count_vec[0,5] += 1
        
    else:
        for iii in string_var:
            if 'China' in iii:
                count_vec[0,0] += 1
            
            if 'U' in iii and 'S' in iii and 'A' in iii:
                count_vec[0,1] += 1
            
            if 'U' in iii and 'K' in iii:
                count_vec[0,2] += 1
            
            if 'Germany' in iii or 'France' in iii or 'Italy' in iii or 'Spain' in iii:
                count_vec[0,3] += 1
            
            if 'Japan' in iii or 'Korea' in iii:
                count_vec[0,4] += 1
            
            if 'India' in iii:
                count_vec[0,5] += 1
            
    return count_vec


# import time
# start = time.process_time()
line = 123; count = 0;
fp = open('data/kg.nt.original', 'r', encoding='utf-8')
newlist = list();
_ = [newlist.append('') for iii in range(N)]
nationlist = np.zeros([N,6], dtype=int)


while line:
    if count % 100000 == 0:
        print(count, end=' ', flush=True);
    
    line = fp.readline()
    x = line.split('> <')
    
    if len(x)<=2:
        continue
    
    temp1 = x[0]+'>'
    temp_act = x[1];
    temp2 = '<'+x[2]
    temp2 = temp2[:-3]
    
    if temp1[0] != '<' or temp2[-1] != '>':
        break
    
    if kg_act[temp_act] == 0 and temp1 in paper_dict:
        line_number = paper_dict[temp1]
        newlist[line_number] = zz_name_dict[temp2]
    
    if kg_act[temp_act] == 2 and temp1 in paper_dict:
        line_number = paper_dict[temp1]
        if temp2 in person_dict:
            author = temp2
            if author in person_to_org_dict:
                org_list = person_to_org_dict[author]
                if type(org_list) is not list:
                    if org_list in city_to_nation_dict:
                        nation_list = city_to_nation_dict[org_list]
                        nationlist[line_number,:] = nationlist[line_number,:] + (Func_nation_count(nation_list)[0,:]>0)
                else:
                    for uuu in org_list:
                        if uuu in city_to_nation_dict:
                            nation_list = city_to_nation_dict[uuu]
                            nationlist[line_number,:] = nationlist[line_number,:] + (Func_nation_count(nation_list)[0,:]>0)
    
    count += 1

fp.close()
# endtime = time.process_time()
# print(endtime - start)

All_Nodes.type = newlist
All_Nodes.China = nationlist[:,0]
All_Nodes.USA   = nationlist[:,1]
All_Nodes.UK    = nationlist[:,2]
All_Nodes.EU    = nationlist[:,3]
All_Nodes.JapanKorea = nationlist[:,4]
All_Nodes.India = nationlist[:,5]


All_Nodes.to_csv('./data/kg_Table.txt')





