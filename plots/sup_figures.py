#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 17 12:15:24 2024

@author: kanep
"""

import os
import csv
import pandas as pd
import numpy as np
import scipy.io
import seaborn as sns
import matplotlib.pyplot as plt

# Sup figure 1
os.chdir('/home/kanep/kg98_scratch/Kane/UCLA/')

file_path = 'sub_list_full.csv'
subjects = []
with open(file_path, 'r') as csvfile:
    reader = csv.reader(csvfile)
    for row in reader:
        # Assuming each row contains only one subject ID
        subject_id = row[0]
        subjects.append(subject_id)

fake_matrix=np.ones((300,300))
ut_mask=np.triu(fake_matrix,k=1)
ut_mask=ut_mask.astype(bool)
all_corrs=np.zeros([9,len(subjects)])
for i,sub in enumerate(subjects):
    matrix_reference = np.array(pd.read_csv('{}/func/AROMA1/FC_Schaef300_AROMA_24P_2P.csv'.format(sub), header=None, index_col=None))
    np.fill_diagonal(matrix_reference,1)
    matrices_to_compare = [np.array(pd.read_csv('{}/func/AROMA{}/FC_Schaef300_AROMA_24P_2P.csv'.format(sub,i), header=None,index_col=None)) for i in range(2, 11)]
    corr_vals=np.zeros(9)
    for k in range(0,9):
        np.fill_diagonal(matrices_to_compare[k],1)
        curr_compare=matrices_to_compare[k]
        
        ref_ut= matrix_reference[ut_mask]
        curr_comp_ut=curr_compare[ut_mask]
        corr_vals[k] = np.corrcoef(ref_ut, curr_comp_ut)[0,1]
    all_corrs[:,i]=corr_vals
            
p=sns.displot(all_corrs.T, kde=True,legend=False)
plt.xlabel('correlation (r)')
plt.legend(labels=["Iteration 2", "Iteration 3", "Iteration 4", "Iteration 5", "Iteration 6", "Iteration 7", 
                   "Iteration 8", "Iteration 9", "Iteration 10"])
plt.title('FC-FC correlations across AROMA Iterations')

colors = sns.color_palette("husl", 9)

# Plotting each distribution in a separate subplot with different colors
fig, axes = plt.subplots(3, 3, figsize=(15, 10))
axes = axes.flatten()

for i in range(9):
    sns.histplot(all_corrs[i, :], kde=False, ax=axes[i], color=colors[i])
    axes[i].set_title(f'Iteration {i + 2}')
    axes[i].set_xlabel('correlation (r)')
    axes[i].set_ylabel('Frequency')
    
    axes[i].spines['top'].set_visible(False)
    axes[i].spines['right'].set_visible(False)

plt.suptitle('FC-FC correlations across AROMA Iterations', fontsize=16)
plt.tight_layout(rect=[0, 0, 1, 0.95])
plt.show()

# Sup Figure 2
# Get performance from LOO results text file in the FIX/Train_data folder 
LOO_res=pd.DataFrame()
LOO_res['TP']=[98.8,98.2,97.1,96.6,94.7,93.2,91.7,89.8]
LOO_res['TN']=[22.7,33.8,53.1,64.3,76.3,85.1,89.0,91.2]
LOO_res['thresh']=[1,2,5,10,20,30,40,50]

sns.set_palette("Pastel2")
plt.figure(figsize=(10, 6))
bar_width = 0.35
r1 = range(len(LOO_res['thresh']))
r2 = [x + bar_width for x in r1]

plt.bar(r1, LOO_res['TP'], width=bar_width, edgecolor='grey', label='TP')
plt.bar(r2, LOO_res['TN'], width=bar_width, edgecolor='grey', label='TN')

plt.xlabel('Threshold')
plt.xticks([r + bar_width/2 for r in range(len(LOO_res['thresh']))], LOO_res['thresh'])
plt.ylabel('Percent (%)')
plt.title('True Positives (TP) and True Negatives (TN) across ICA-FIX classifier thresholds')
plt.legend()
plt.grid(axis='y', linestyle='--', alpha=0.7)
plt.show()

# Sup Figure 3
HCP_opt_k=np.array(pd.read_csv('/home/kanep/kg98_scratch/Kane/HCP/Results/Modularity/Optimal_iterations_consensus_2P.csv', header=None))
GSP_FIX_opt_k=np.array(pd.read_csv('/home/kanep/kg98_scratch/Kane/GSP/FIX/Results/Modularity/Optimal_iterations_consensus_2P.csv', header=None))
GSP_AROMA_opt_k=np.zeros([10,1568])
UCLA_opt_k=np.zeros([10,121])
for i,k in enumerate(range(1,11)):
    GSP_tmp=np.array(pd.read_csv('/home/kanep/kg98_scratch/Kane/GSP/Results/Modularity/Optimal_iterations_consensus_2P_AROMA{}.csv'.format(k), header=None))
    UCLA_tmp=np.array(pd.read_csv('/home/kanep/kg98_scratch/Kane/UCLA/Results/Modularity/Optimal_iterations_consensus_2P_AROMA{}.csv'.format(k), header=None))
    GSP_AROMA_opt_k[i,:]=GSP_tmp[0,:]
    UCLA_opt_k[i,:]=UCLA_tmp[0,:]
    
plt.figure()
sns.set_palette('Accent')
sns.violinplot(HCP_opt_k.T, legend=False)
plt.title('Optimal DiCER Iterations (HCP)')
plt.ylabel('DiCER Iteration')
plt.xticks('HCP')
sns.despine(top=True, right=True)

plt.figure()
sns.set_palette('Accent')
sns.violinplot(GSP_FIX_opt_k.T, legend=False)
plt.title('Optimal DiCER Iterations (GSP FIX)')
plt.ylabel('DiCER Iteration')
plt.xticks('GSP')
sns.despine(top=True, right=True)


plt.figure()
sns.violinplot(UCLA_opt_k.T, legend=False)
plt.title('Optimal DiCER Iterations (CNP)')
plt.ylabel('DiCER Iteration')
plt.xticks(ticks=range(10), labels=["1", "2", "3", "4", "5", 
                                     "6", "7", "8", "9", "10"])
plt.xlabel("AROMA Iteration")
plt.ylim([0,5.5])
sns.despine(top=True, right=True)

plt.figure()
sns.violinplot(GSP_AROMA_opt_k.T, legend=False)
plt.title('Optimal DiCER Iterations (GSP AROMA)')
plt.ylabel('DiCER Iteration')
plt.xticks(ticks=range(10), labels=["1", "2", "3", "4", "5", 
                                     "6", "7", "8", "9", "10"])
plt.xlabel("AROMA Iteration")
plt.ylim([0,5.5])
sns.despine(top=True, right=True)


# Sup Figure 4
os.chdir('/home/kanep/kg98_scratch/Kane/UCLA')
QC_corrs=scipy.io.loadmat('Results/QC/QC_FC_corrs.mat')
QC_corrs=QC_corrs['collated_QC']
preprocs=np.array(pd.read_csv("Results/QC/preprocs.csv", header=None))
a_iters=10
n_edges=QC_corrs['AROMA1']; n_edges=n_edges[0,0]; n_edges=n_edges.shape[1]

# Get mean QC-FC Distributions
curr_proc=np.zeros([a_iters,n_edges])
mean_qc=np.zeros([len(preprocs),n_edges])
min_qc=np.zeros([len(preprocs),n_edges])
max_qc=np.zeros([len(preprocs),n_edges])
for p in range(len(preprocs)):
    for a in range(a_iters):
        curr=QC_corrs['AROMA{}'.format(a+1)]
        curr=curr[0,0]
        curr_proc[a,:]=curr[p,:]
    mean_qc[p,:]=np.mean(curr_proc, axis=0)
    min_qc[p,:]=np.min(curr_proc, axis=0)
    max_qc[p,:]=np.max(curr_proc, axis=0)

mmmm_df=pd.DataFrame()
mmmm_df['min']=np.mean(min_qc,axis=1)
mmmm_df['average']=np.mean(mean_qc,axis=1)
mmmm_df['max']=np.mean(max_qc,axis=1)
mmmm_df['preprocs']=preprocs[:,0]

sns.scatterplot(mmmm_df)
plt.xticks(ticks=(range(len(mmmm_df))),labels=mmmm_df['preprocs'], rotation=75)
plt.ylabel('QC-FC (r)')
plt.title('Differences in QC-FC across AROMA Iterations (GSP)')
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)

# Sup Figure 5
# For FIX
os.chdir('/home/kanep/kg98_scratch/Kane/GSP/FIX')
sns.set_palette('Accent')
file_path = '../sub_list.txt'
subjects = []
with open(file_path, 'r') as csvfile:
    reader = csv.reader(csvfile)
    for row in reader:
        # Assuming each row contains only one subject ID
        subject_id = row[0]
        subjects.append(subject_id)

n_regs=10
all_corrs=np.zeros([n_regs,len(subjects)])
missing_subs=[]
for i,sub in enumerate(subjects):
    if os.path.exists('{}/{}_dbscan_liberal_regressors.csv'.format(sub,sub)):
        dic_regs=np.array(pd.read_csv('{}/{}_dbscan_liberal_regressors.csv'.format(sub,sub), index_col=0))
        GS=np.array(pd.read_csv('{}/GS_ts.txt'.format(sub), header=None))
        curr_corrs=np.zeros(len(dic_regs))
        for r in range(len(dic_regs)):
            curr_reg=dic_regs[r,:]
            curr_corrs[r]=np.corrcoef(curr_reg,GS[:,0])[0,1]
        all_corrs[:,i]=curr_corrs
    else:
        missing_subs.append(i)
    print(sub)
all_corrs=np.delete(all_corrs, missing_subs, axis=1)
    
sns.violinplot(all_corrs.T)
plt.xticks(range(0, 10), labels=range(1, 11))
plt.xlabel('DiCER Iteration')
plt.ylabel('Correlation with Global Signal (r)')
plt.title('GSP FIX')

# For AROMA
os.chdir('/home/kanep/kg98_scratch/Kane/GSP')
n_regs=10
a_iters=10
all_corrs=np.zeros([a_iters,n_regs,len(subjects)])
missing_subs=[]
for a in range(1,a_iters+1):
    for i,sub in enumerate(subjects):
        if os.path.exists('{}/func/AROMA{}/{}_dbscan_liberal_regressors.csv'.format(sub,a,sub)):
            dic_regs=np.array(pd.read_csv('{}/func/AROMA{}/{}_dbscan_liberal_regressors.csv'.format(sub,a,sub), index_col=0))
            GS=np.array(pd.read_csv('{}/func/AROMA{}/GS_ts.txt'.format(sub,a), header=None))
            curr_corrs=np.zeros(len(dic_regs))
            for r in range(len(dic_regs)):
                curr_reg=dic_regs[r,:]
                curr_corrs[r]=np.corrcoef(curr_reg,GS[:,0])[0,1]
                all_corrs[a-1,:,i]=curr_corrs
        else:
            missing_subs.append(i)
        print(sub)
all_corrs=np.delete(all_corrs, missing_subs, axis=2)
mean_corrs = np.mean(all_corrs, axis=0)

sns.violinplot(mean_corrs.T)
plt.xticks(range(0, 10), labels=range(1, 11))
plt.xlabel('DiCER Iteration')
plt.ylabel('Correlation with Global Signal (r)')
plt.title('GSP AROMA')