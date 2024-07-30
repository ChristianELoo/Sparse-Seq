#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created: 09/06/2023 
Author: Jojo Zhu

Description: 

"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from operator import add

data = pd.read_csv('downsampling_data_ALL.csv')
full_stats = pd.read_csv('full_bam_statistics.csv')

# create dictionary with key = GSM number and value = perc total modification from experiment full stats
full_stats_dic = {}
for x in range(0,len(full_stats)):
    sample = full_stats['gsm'][x]
    mod = full_stats['modC_perc'][x]
    full_stats_dic[f'{sample}']=mod
    
# create list with total modification from full dataset for each row in downsampled data set according to the experiment
# using dictionary created above
# add that list as a row to the data dataframe
full_stats_list = []
for gsm in data['gsm']:
    full_stats_list.append(full_stats_dic[f'{gsm}'])
    
data['full_stat']=full_stats_list

data['error']=data['full_stat']-data['mC_perc']

# create dictionary w/ each experiment as key and its corresponsding data frame as value
bs = data.groupby('exp').get_group('BS')
tab = data.groupby('exp').get_group('TAB')
ace = data.groupby('exp').get_group("ACE")
em_mouse = data.groupby(['exp','cell']).get_group(("EM",'Mouse neurons'))
em_human = data.groupby(['exp','cell']).get_group(("EM",'human_bcells'))

experiments = {'bs':bs,'tab':tab,'ace':ace,'em_mouse':em_mouse,'em_human':em_human}

# coverages list
coverages = {}

for key in experiments:
    coverage = []
    for name, group in experiments[key].groupby('coverage'):
        coverage.append(name)
    coverages[key]=coverage
    
#%%
# create group by experiment
groups_exp_BS = experiments['bs'].groupby('gsm')
groups_exp_TAB = experiments['tab'].groupby('gsm')

# iterate over experiment groups...
experiments_BS = {}
experiments_TAB = {}

# starting with BS data
for name_exp, group in groups_exp_BS:
    #...and further group by coverage PERCENT
    group_coverage_perc = group.groupby('coverage_perc')
    # create dictionary for each coverage percent group's mC data
    boxes = {}
    # iterate over coverage percent groups...
    for name_coverage_perc, coverage_perc_group in group_coverage_perc:
        #... create list from mC_perc
        x = coverage_perc_group['mC_perc'].tolist()
        # and save that list, under its coverage percent in the boxes dictionary
        boxes[f'{name_coverage_perc}']=x
    # save the boxes dictionary in the experiments dictionay under its experiment name
    experiments_BS[f'{name_exp}']=boxes
    
# repeat for TAB data
for name_exp, group in groups_exp_TAB:
    group_coverage_perc = group.groupby('coverage_perc')
    boxes = {}
    for name_coverage_perc, coverage_perc_group in group_coverage_perc:
        x = coverage_perc_group['mC_perc'].tolist()
        boxes[f'{name_coverage_perc}']=x
    experiments_TAB[f'{name_exp}']=boxes
    
# create dictionary with just full stats from BS and TAB experiments
BS_TAB_stats = {}

for key in experiments_BS:
    BS_TAB_stats[f'{key}']=[full_stats_dic[f'{key}']]

for key in experiments_TAB:
    BS_TAB_stats[f'{key}']=[full_stats_dic[f'{key}']]

#%%%
plt.style.use('ggplot')
plt.rcParams.update({'font.family':'Arial',"axes.spines.bottom":True})

fig, ax = plt.subplots(figsize=(9,3.7))

ax.set_facecolor('white')
ax.tick_params(grid_color='gray',grid_alpha=0)

colors = {'GSM1180315':'tab:red','GSM1180316':"tab:red","GSM1180317":"tab:red",
         "GSM1541958":"tab:orange","GSM1541959":"tab:orange",
         "GSM1173794":"tab:green",
         "GSM1173795":"tab:blue",
         "GSM1180306":"tab:purple","GSM1180307":"tab:purple",'GSM1180308':"tab:purple"}

position_list = [1,4,7,10,13,16,19,22,25,28,31,34,37]
counter = -0.8

# add horizontal line for each experiment corresponding to the modified cytosine level from the full (nondownsampled) experiment
for key in BS_TAB_stats:
    ax.axhline(y=BS_TAB_stats[f'{key}'][0],color=colors[f'{key}'],ls='--',lw=0.7)

# add boxplots for downsampling data at each coverage level for BS-seq experiments
for key in experiments_BS:
    boxes = experiments_BS[f'{key}']
    labels, points = [*zip(*boxes.items())]
    new_pos = [z+counter for z in position_list]
    bp = ax.boxplot(points,widths=0.4,positions=new_pos,patch_artist=True,showcaps=False,showfliers=False,
                   boxprops=dict(facecolor=colors[f'{key}'],color='black'),
                   capprops=dict(color='black'),
                   whiskerprops=dict(color='black'),
                   flierprops=dict(color='black', markeredgecolor='black'),
                   medianprops=dict(color='black'),manage_ticks=False
                   )
    counter+=0.4

position_list = [1,4,7,10,13,16,19,22,25,28,31,34,37]
counter = -0.8

# add boxplots for downsampling data at each coverage level for TAB-seq experiments
for key in experiments_TAB:
    boxes = experiments_TAB[f'{key}']
    labels, points = [*zip(*boxes.items())]
    new_pos = [z+counter for z in position_list]
    bp = ax.boxplot(points,widths=0.4,positions=new_pos,patch_artist=True,showcaps=False,showfliers=False,
                   boxprops=dict(facecolor=colors[f'{key}'],color='black'),
                   capprops=dict(color='black'),
                   whiskerprops=dict(color='black'),
                   flierprops=dict(color='black', markeredgecolor='black'),
                   medianprops=dict(color='black'),manage_ticks=False
                   )
    counter+=0.4

# set x ticks to appropriate coverage levels
tick_list = [1,4,7,10,13,16,19,22,25,28,31,34,37]    
ax.set_xticks(ticks=tick_list)
ax.set_xticklabels([f'{i*100:.5f}' for i in coverages['bs']],fontsize=8)

ax.set_yticks(np.arange(0,9,1))

plt.tick_params(axis='y',which='major',labelsize=8,color='black',labelcolor='black')
plt.tick_params(axis='x',which='major',color='black',labelcolor='black')

ax.set_ylabel('total cytosine modification (%)',fontsize=12,c='black')
ax.set_xlabel('genomic coverage (%)',fontsize=12,c='black')
ax.spines['bottom'].set_color('black')
ax.spines['left'].set_color('black')


plt.show()

#%%
# create empty dictionary to store data
stats_BS_TAB = {}

# for each experiment processed...
for name, group in groups_exp_BS:
    #...separated data into group by level (coverage) they were downsampled to
    cov_groups = group.groupby('coverage')
    cov_dict = {}
    # for each coverage (downsampled) level group of replicates...
    for cov, groupC in cov_groups:
        #...calculate the average value for modified cytosine...
        avg = groupC['mC_perc'].mean()
        #...calculate the standard deviation for modified cytosine...
        std = groupC['mC_perc'].std()
        #...calculate the bias...
        bias = abs(groupC['error'].mean())
        #...calculate the coefficient of variation...
        cov = std/avg
        #...calculate the total analytical error...
        tae = bias*100 + 1.96*100*cov
        #...pull out the coverage (this "mean" is just a means to pull the number out--each downsampled replicate will have the same number for coverage percent)...
        covg_perc = groupC['coverage_perc'].mean()
        #...pull out the modified C number from the full, nondownsampled data set (this number is also the same for each downsampled replicate)...
        full_perc = groupC['full_stat'].mean()
        #...compile statistics calculated above into a list...
        stat_list = [std,cov,avg,bias,tae,covg_perc,full_perc]
        #...add list of statistics to coverage dictionary under key corresponding to the coverage value to which the replicates were downsampled
        cov_dict[f'{covg_perc}']=stat_list
    # add dictionary of statistics list (value) by coverage level (key) as a value to dictionary where the key is the experiment (GSM)
    stats_BS_TAB[f'{name}']=cov_dict

# repeat as above but for TAB-seq data
for name, group in groups_exp_TAB:
    cov_groups = group.groupby('coverage')
    cov_dict = {}
    for cov, groupC in cov_groups:
        avg = groupC['mC_perc'].mean()
        std = groupC['mC_perc'].std()
        bias = abs(groupC['error'].mean())
        cov = std/avg
        tae = bias*100 + 1.96*100*cov
        covg_perc = groupC['coverage_perc'].mean()
        full_perc = groupC['full_stat'].mean()
        stat_list = [std,cov,avg,bias,tae,covg_perc,full_perc]
        cov_dict[f'{covg_perc}'] = stat_list
    stats_BS_TAB[f'{name}'] = cov_dict

# indeces in each statistics dictionary (stats_BS_TAB['GSM']['coverage level][i]) correspond to:
# 0 = standard deviation
# 1 = coefficient of variance
# 2 = average
# 3 = bias
# 4 = total analytical error
# 5 = % coverage
# 6 = experiment full dataset average

#%%
plt.style.use('ggplot')
plt.rcParams.update({'font.family':'Arial'})

fig, ax = plt.subplots(figsize=(8,4))

ax.set_facecolor('white')
ax.tick_params(grid_color='gray',grid_alpha=0.3)

colors = {'GSM1180315':'tab:red','GSM1180316':"tab:red","GSM1180317":"tab:red",
         "GSM1541958":"tab:orange","GSM1541959":"tab:orange",
         "GSM1173794":"tab:green",
         "GSM1173795":"tab:blue",
         "GSM1180306":"tab:purple","GSM1180307":"tab:purple",'GSM1180308':"tab:purple"}

for key in stats_BS_TAB:
    for power in stats_BS_TAB[f'{key}']:
        ax.plot(stats_BS_TAB[f'{key}'][f'{power}'][-2],stats_BS_TAB[f'{key}'][f'{power}'][4],ls='',marker='.',color=colors[f'{key}'],markersize=7)
        
ax.set_xscale('log')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

tick_list = [i*100 for i in coverages['bs']] 
ax.set_xticks(ticks=tick_list)
ax.set_xticklabels([f'{i*100:.5f}' for i in coverages['bs']],fontsize=8)

#ax.axhline(y=2)
plt.tick_params(axis='x',which='major',labelsize=7)
plt.tick_params(axis='y',which='major',labelsize=7)

plt.savefig('240105_BS_TAB_TAE.pdf',bbox_inches='tight',pad_inches=0)


plt.show()

#%%
# export data used in above plot to CSV

bs_tab_tae_data = pd.DataFrame()

gsm = []
power_list = []
coverage = []
tae = []
full_dataset_mod = []

# for each experiment (key)...
for key in stats_BS_TAB:
    #...iterate over coverage levels in downsampling...
    for power in stats_BS_TAB[f'{key}']:
        # add experiment to GSM list
        gsm.append(key)
        # add power of downsampling to power list
        power_list.append(power)
        # add coverage of downsampling to coverage list
        coverage.append(stats_BS_TAB[f'{key}'][f'{power}'][-2])
        # add TAE to tae list
        tae.append(stats_BS_TAB[f'{key}'][f'{power}'][4])
        #add full dataset (nondownsampled) modified cytosine information to full_dataset_mod list
        full_dataset_mod.append(stats_BS_TAB[f'{key}'][f'{power}'][6])
        
# add lists created above to dataframe
bs_tab_tae_data["gsm"]=gsm
#bs_tab_tae_data['power']=power_list
bs_tab_tae_data['coverage_perc']=coverage
bs_tab_tae_data['tae']=tae
bs_tab_tae_data['full_dataset_mod']=full_dataset_mod

# add columns for experiment pipeline and cell type
exp=[]
cell=[]
for i in range(0,len(bs_tab_tae_data)):
    #print(data.iloc[i]['gsm'])
    if '118031' in bs_tab_tae_data.iloc[i]['gsm']:
        exp.append('BS')
        cell.append('Mouse ESCs')
    if '154' in bs_tab_tae_data.iloc[i]['gsm']:
        exp.append('BS')
        cell.append('Mouse neurons')
    if '118030' in bs_tab_tae_data.iloc[i]['gsm']:
        exp.append('TAB')
        cell.append("Mouse ESCs")
    if '320718' in bs_tab_tae_data.iloc[i]['gsm']:
        exp.append('ACE')
    if '3207181' in bs_tab_tae_data.iloc[i]['gsm']:
        cell.append('Mouse ESCs')
    if '3207182' in bs_tab_tae_data.iloc[i]['gsm']:
        cell.append('Mouse ESCs')
    if '3207183' in bs_tab_tae_data.iloc[i]['gsm']:
        cell.append('TetTKO mESCs')
    if '3207184' in bs_tab_tae_data.iloc[i]['gsm']:
        cell.append('TetTKO mESCs')
    if '3207185' in bs_tab_tae_data.iloc[i]['gsm']:
        cell.append('Mouse neurons')
    if '3207186' in bs_tab_tae_data.iloc[i]['gsm']:
        cell.append('Mouse neurons')
    if '1173794' in bs_tab_tae_data.iloc[i]['gsm']:
        cell.append('Mouse fetal cortex')
        exp.append('TAB')
    if '1173795' in bs_tab_tae_data.iloc[i]['gsm']:
        cell.append('Mouse adult cortex')
        exp.append('TAB')
    if 'NA12878' in bs_tab_tae_data.iloc[i]['gsm']:
        cell.append('human_bcells')
        exp.append('EM')
    if "584964" in bs_tab_tae_data.iloc[i]['gsm']:
        cell.append('Mouse neurons')
        exp.append('EM')
        
bs_tab_tae_data['exp']=exp
bs_tab_tae_data['cell']=cell

# create group by experiment
groups_exp_ACE = experiments['ace'].groupby('gsm')
groups_exp_EM_m = experiments['em_mouse'].groupby('gsm')
groups_exp_EM_h = experiments['em_human'].groupby('gsm')

# iterate over experiment groups...
experiments_ACE = {}
experiments_EM_m = {}
experiments_EM_h = {}

# starting with ACE data
for name_exp, group in groups_exp_ACE:
    #...and further group by coverage PERCENT
    group_coverage_perc = group.groupby('coverage_perc')
    # create dictionary for each coverage percent group's mC data
    boxes = {}
    # iterate over coverage percent groups...
    for name_coverage_perc, coverage_perc_group in group_coverage_perc:
        #... create list from mC_perc
        x = coverage_perc_group['mC_perc'].tolist()
        # and save that list, under its coverage percent in the boxes dictionary
        boxes[f'{name_coverage_perc}']=x
    # save the boxes dictionary in the experiments dictionay under its experiment name
    experiments_ACE[f'{name_exp}']=boxes

# repeat for EM mouse data
for name_exp, group in groups_exp_EM_m:
    group_coverage_perc = group.groupby('coverage_perc')
    boxes = {}
    for name_coverage_perc, coverage_perc_group in group_coverage_perc:
        x = coverage_perc_group['mC_perc'].tolist()
        boxes[f'{name_coverage_perc}']=x
    experiments_EM_m[f'{name_exp}']=boxes

# repeat again for EM human data
for name_exp, group in groups_exp_EM_h:
    group_coverage_perc = group.groupby('coverage_perc')
    boxes = {}
    for name_coverage_perc, coverage_perc_group in group_coverage_perc:
        x = coverage_perc_group['mC_perc'].tolist()
        boxes[f'{name_coverage_perc}']=x
    experiments_EM_h[f'{name_exp}']=boxes
    
# create dictionary with just full stats from BS and TAB experiments
ACE_EM_stats = {}

for key in experiments_ACE:
    ACE_EM_stats[f'{key}']=[full_stats_dic[f'{key}']]

for key in experiments_EM_m:
    ACE_EM_stats[f'{key}']=[full_stats_dic[f'{key}']]

for key in experiments_EM_h:
    ACE_EM_stats[f'{key}']=[full_stats_dic[f'{key}']]
    
plt.style.use('ggplot')
plt.rcParams.update({'font.family':'Arial',"axes.spines.bottom":True,"axes.spines.top":True})
fig, ax = plt.subplots(figsize=(9,3.7))

ax.set_facecolor('white')
ax.tick_params(grid_color='gray',grid_alpha=0)

colors = {"GSM5849646":'tab:brown','GSM5849647':'tab:brown',
        'NA12878_50ng_rep1':"tab:gray",'NA12878_50ng_rep2':'tab:gray',
        'GSM3207181':'tab:pink','GSM3207182':"tab:pink",
        'GSM3207183':'tab:olive','GSM3207184':'tab:olive',
        'GSM3207185':'tab:cyan','GSM3207186':'tab:cyan'}

position_list = [1,4,7,10,13,16,19,22,25,28,31,34,37]
counter = -0.8

for key in experiments_EM_m:
    boxes = experiments_EM_m[f'{key}']
    labels, points = [*zip(*boxes.items())]
    new_pos = [z+counter for z in position_list]
    bp = ax.boxplot(points,widths=0.4,positions=new_pos,patch_artist=True,showcaps=False,showfliers=False,
                   boxprops=dict(facecolor=colors[f'{key}'],color='black'),
                   capprops=dict(color='black'),
                   whiskerprops=dict(color='black'),
                   flierprops=dict(color='black',markeredgecolor='black'),
                   medianprops=dict(color='black'),manage_ticks=False)
    counter+=0.4

for key in experiments_EM_h:
    boxes = experiments_EM_h[f'{key}']
    labels, points = [*zip(*boxes.items())]
    new_pos = [z+counter for z in position_list]
    bp = ax.boxplot(points,widths=0.4,positions=new_pos,patch_artist=True,showcaps=False,showfliers=False,
                   boxprops=dict(facecolor=colors[f'{key}'],color='black'),
                   capprops=dict(color='black'),
                   whiskerprops=dict(color='black'),
                   flierprops=dict(color='black',markeredgecolor='black'),
                   medianprops=dict(color='black'),manage_ticks=False)
    counter+=0.4
    
position_list = [1,4,7,10,13,16,19,22,25,28,31,34,37]
counter = -0.8

for key in experiments_ACE:
    boxes = experiments_ACE[f'{key}']
    labels, points = [*zip(*boxes.items())]
    new_pos = [z+counter for z in position_list]
    bp = ax.boxplot(points,widths=0.4,positions=new_pos,patch_artist=True,showcaps=False,showfliers=False,
                   boxprops=dict(facecolor=colors[f'{key}'],color='black'),
                   capprops=dict(color='black'),
                   whiskerprops=dict(color='black'),
                   flierprops=dict(color='black',markeredgecolor='black'),
                   medianprops=dict(color='black'),manage_ticks=False)
    counter+=0.4
    
for key in ACE_EM_stats:
    ax.axhline(y=ACE_EM_stats[f'{key}'][0],color=colors[f'{key}'],ls='--',lw=0.7)
    
tick_list = [1,4,7,10,13,16,19,22,25,28,31,34,37]    
ax.set_xticks(ticks=tick_list)
# labels for ACE-seq
ax.set_xticklabels([f'{i*100:.5f}' for i in coverages['ace']],
                  fontsize=8)

secax = ax.secondary_xaxis('top')
#secax.set_xlabel('genomic coverage (%) EM-Seq human',fontsize=12,c='black')
secax.set_xticks(ticks=tick_list)
secax.set_xticklabels([f'{i*100:.5f}' for i in coverages['em_human']],fontsize=8)
secax.tick_params(axis='x',which='major',labelsize=8,color='black',labelcolor='black')


ax.set_yticks(np.arange(0,9,1))

#ax.set_ylabel('total cytosine modification (%)',fontsize=12,c='black')
#ax.set_xlabel('genomic coverage (%) ACE-Seq/EM-Seq mouse',fontsize=12,c='black')
ax.spines['bottom'].set_color('black')
ax.spines['left'].set_color('black')
ax.spines['top'].set_color('black')


plt.tick_params(axis='y',which='major',labelsize=8,color='black',labelcolor='black')
plt.tick_params(axis='x',which='major',color='black',labelcolor='black')

plt.show()

#%% Separated EM Plot 
plt.style.use('ggplot')
plt.rcParams.update({'font.family':'Arial',"axes.spines.bottom":True,"axes.spines.top":True})
fig, ax = plt.subplots(figsize=(9,3.7))

ax.set_facecolor('white')
ax.tick_params(grid_color='gray',grid_alpha=0)

colors = {"GSM5849646":'tab:purple','GSM5849647':'tab:purple',
        'NA12878_50ng_rep1':"tab:blue",'NA12878_50ng_rep2':'tab:blue',
        'GSM3207181':'tab:pink','GSM3207182':"tab:pink",
        'GSM3207183':'tab:olive','GSM3207184':'tab:olive',
        'GSM3207185':'tab:cyan','GSM3207186':'tab:cyan'}

position_list = [1,4,7,10,13,16,19,22,25,28,31,34,37]
counter = -0.8

for key in experiments_EM_m:
    boxes = experiments_EM_m[f'{key}']
    labels, points = [*zip(*boxes.items())]
    new_pos = [z+counter for z in position_list]
    bp = ax.boxplot(points,widths=0.4,positions=new_pos,patch_artist=True,showcaps=False,showfliers=False,
                   boxprops=dict(facecolor=colors[f'{key}'],color='black'),
                   capprops=dict(color='black'),
                   whiskerprops=dict(color='black'),
                   flierprops=dict(color='black',markeredgecolor='black'),
                   medianprops=dict(color='black'),manage_ticks=False)
    counter+=0.4
    ax.axhline(y=ACE_EM_stats[f'{key}'][0],color=colors[f'{key}'],ls='--',lw=0.7)

for key in experiments_EM_h:
    boxes = experiments_EM_h[f'{key}']
    labels, points = [*zip(*boxes.items())]
    new_pos = [z+counter for z in position_list]
    bp = ax.boxplot(points,widths=0.4,positions=new_pos,patch_artist=True,showcaps=False,showfliers=False,
                   boxprops=dict(facecolor=colors[f'{key}'],color='black'),
                   capprops=dict(color='black'),
                   whiskerprops=dict(color='black'),
                   flierprops=dict(color='black',markeredgecolor='black'),
                   medianprops=dict(color='black'),manage_ticks=False)
    counter+=0.4
    ax.axhline(y=ACE_EM_stats[f'{key}'][0],color=colors[f'{key}'],ls='--',lw=0.7)
    
position_list = [1,4,7,10,13,16,19,22,25,28,31,34,37]
counter = -0.8

tick_list = [1,4,7,10,13,16,19,22,25,28,31,34,37]    
ax.set_xticks(ticks=tick_list)
# labels for ACE-seq
ax.set_xticklabels([f'{i*100:.5f}' for i in coverages['em_human']],
                  fontsize=8)
secax.tick_params(axis='x',which='major',labelsize=8,color='black',labelcolor='black')


ax.set_yticks(np.arange(0,9,1))

#ax.set_ylabel('total cytosine modification (%)',fontsize=12,c='black')
#ax.set_xlabel('genomic coverage (%) ACE-Seq/EM-Seq mouse',fontsize=12,c='black')
ax.spines['bottom'].set_color('black')
ax.spines['left'].set_color('black')
ax.spines['top'].set_color('black')


plt.tick_params(axis='y',which='major',labelsize=8,color='black',labelcolor='black')
plt.tick_params(axis='x',which='major',color='black',labelcolor='black')

plt.show()

#%% Separated ACE plot 
plt.style.use('ggplot')
plt.rcParams.update({'font.family':'Arial',"axes.spines.bottom":True,"axes.spines.top":True})
fig, ax = plt.subplots(figsize=(9,3.7))

ax.set_facecolor('white')
ax.tick_params(grid_color='gray',grid_alpha=0)

colors = {"GSM5849646":'tab:purple','GSM5849647':'tab:purple',
        'NA12878_50ng_rep1':"tab:blue",'NA12878_50ng_rep2':'tab:blue',
        'GSM3207181':'tab:pink','GSM3207182':"tab:pink",
        'GSM3207183':'tab:olive','GSM3207184':'tab:olive',
        'GSM3207185':'tab:cyan','GSM3207186':'tab:cyan'}

position_list = [1,4,7,10,13,16,19,22,25,28,31,34,37]
counter = -0.8

for key in experiments_ACE:
    boxes = experiments_ACE[f'{key}']
    labels, points = [*zip(*boxes.items())]
    new_pos = [z+counter for z in position_list]
    bp = ax.boxplot(points,widths=0.4,positions=new_pos,patch_artist=True,showcaps=False,showfliers=False,
                   boxprops=dict(facecolor=colors[f'{key}'],color='black'),
                   capprops=dict(color='black'),
                   whiskerprops=dict(color='black'),
                   flierprops=dict(color='black',markeredgecolor='black'),
                   medianprops=dict(color='black'),manage_ticks=False)
    counter+=0.4
    ax.axhline(y=ACE_EM_stats[f'{key}'][0],color=colors[f'{key}'],ls='--',lw=0.7)
    
    
tick_list = [1,4,7,10,13,16,19,22,25,28,31,34,37]    
ax.set_xticks(ticks=tick_list)
# labels for ACE-seq
ax.set_xticklabels([f'{i*100:.5f}' for i in coverages['ace']],
                  fontsize=8)

ax.set_yticks(np.arange(0,9,1))

ax.set_ylim(0, 3)

#ax.set_ylabel('total cytosine modification (%)',fontsize=12,c='black')
#ax.set_xlabel('genomic coverage (%) ACE-Seq/EM-Seq mouse',fontsize=12,c='black')
ax.spines['bottom'].set_color('black')
ax.spines['left'].set_color('black')
ax.spines['top'].set_color('black')


plt.tick_params(axis='y',which='major',labelsize=8,color='black',labelcolor='black')
plt.tick_params(axis='x',which='major',color='black',labelcolor='black')

plt.show()
