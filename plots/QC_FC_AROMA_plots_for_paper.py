#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 27 10:42:41 2023

@author: kanep
"""

import os
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import scipy.io


os.chdir("/home/kanep/kg98_scratch/Kane/UCLA")

colors = [(1, 1, 0.3), (1, 0.95, 0.5), (0.176, 0.839, 0.341), (0.38, 0.969, 0.529), 
          (0.817, 0.188, 0.812), (0.961, 0.431, 0.914), (0.153, 0.729, 0.812), (0.396, 0.875, 0.941), 
          (0.902, 0.576, 0.094), (0.961, 0.773, 0.486), (0.878, 0.071, 0.145), (0.969, 0.38, 0.435), 
          (0.298, 0.31, 0.969), (0.392, 0.4, 0.961),(0.678, 0.659, 0.678), (0.839, 0.824, 0.839), 
          (0.96, 0.15, 0.52), (1, 0.41, 0.7)]
palette=sns.color_palette(colors)

labels = np.array([("1"), ("2"), ("3*"), ("4*"), ("5*"), ("6*"), ("7*"), ("8*"), ("9*"), ("10*"),
          ("11*"), ("12*"), ("13*"), ("14*"), ("15*"), ("16*"), ("17*"), ("18*")])
# Load results
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
    
    mean_vector = np.mean(curr_proc, axis=1)
    mean_index = np.abs(mean_vector - np.mean(mean_vector)).argmin()
    mean_qc[p,:]=curr_proc[mean_index,:]
    min_index=np.argmin(mean_vector)
    min_qc[p,:]=curr_proc[min_index,:]
    max_index=np.argmax(mean_vector)
    max_qc[p,:]=curr_proc[max_index,:]
    
    mean_qc[p,:]=np.mean(curr_proc, axis=0)
    min_qc[p,:]=np.min(curr_proc, axis=0)
    max_qc[p,:]=np.max(curr_proc, axis=0)

# Combine Results
mean_qc=pd.DataFrame(mean_qc)
mean_qc.index=labels
mean_qc=pd.DataFrame(mean_qc.stack())
mean_qc=mean_qc.rename(columns={0: 'mean'})
mean_qc['preprocs']=mean_qc.index.get_level_values(0)

min_qc=pd.DataFrame(min_qc)
min_qc=pd.DataFrame(min_qc.stack())
min_qc=min_qc.rename(columns={0: 'min'})
max_qc=pd.DataFrame(max_qc)
max_qc=pd.DataFrame(max_qc.stack())
max_qc=max_qc.rename(columns={0: 'max'})

mmm_qc=pd.DataFrame(mean_qc)
mmm_qc['min']=np.array(min_qc['min'])
mmm_qc['max']=np.array(max_qc['max'])

# Plot the correlation distributions
sns.set_theme(style="white", rc={"axes.facecolor": (0, 0, 0, 0), 'axes.linewidth':2})
g = sns.FacetGrid(mmm_qc, palette=sns.color_palette(colors), row='preprocs', hue='preprocs', aspect=6, height=2.2)
g.map_dataframe(sns.kdeplot, x='mean', fill=True, alpha=1)
g.map_dataframe(sns.kdeplot, x='mean', color='black')
g.map_dataframe(sns.kdeplot, x='min', fill=False, alpha=0.7)
g.map_dataframe(sns.kdeplot, x='max', fill=False, alpha=0.7)
plt.tight_layout()
def label(x, color, label):
    ax = plt.gca()
    ax.text(0, .2, label, color='black', fontsize=13,
            ha="left", va="center", transform=ax.transAxes)
    ax.axvline(x=0, color='black', linestyle='--')
g.map(label, "preprocs")
g.set_titles("")
g.set(yticks=[])
g.set(ylabel="")
g.set(xlabel="QC-FC (r)")
g.set(xlim=[-0.4, 0.5])
g.despine(left=True)
g.fig.subplots_adjust(hspace=-.5)
#plt.suptitle("", fontsize=20, x=0.2, y=0.95)
plt.subplots_adjust(bottom=0.15)



# Load Distance Dependence
# This block has to be run with seaborn version 0.11 (module load anaconda 2019 version)
# The xerr doesn't play nicely with later versions.
dist_dep=scipy.io.loadmat('Results/QC/distance_dependence.mat')
dist_dep=dist_dep['collated_dist_dep']

# Get mean and min/max distance dependence.
curr_dist=np.zeros(a_iters)
mean_dist=np.zeros(len(preprocs))
min_dist=np.zeros(len(preprocs))
max_dist=np.zeros(len(preprocs))
for p in range(len(preprocs)):
        for a in range(10):
            curr=dist_dep['AROMA{}'.format(a+1)]
            curr=curr[0,0]
            curr_dist[a]=curr[:,p]
        mean_dist[p]=np.mean(curr_dist, axis=0)
        min_dist[p]=np.min(curr_dist, axis=0)
        max_dist[p]=np.max(curr_dist, axis=0)
mmm_distdep=pd.DataFrame()
mmm_distdep['preprocs']=labels
mmm_distdep['min']=min_dist
mmm_distdep['mean']=mean_dist
mmm_distdep['max']=max_dist
xerr=[np.subtract(mmm_distdep['mean'], mmm_distdep['min']), np.subtract(mmm_distdep['max'], mmm_distdep['mean'])]

plt.figure()
ax=sns.barplot(data=mmm_distdep, x='mean', y='preprocs', xerr=xerr, palette=palette)
#ax.bar_label(ax.containers[0], fmt='%.3f')
plt.xlim(-0.3,0.2)
plt.axvline(x=0, color='grey', linestyle='-')
plt.xlabel('Spearmans rho (r)')
plt.ylabel=('')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.subplots_adjust(left=0.4)
#plt.title('QC-FC Distant Dependence')

# Draw VE by PC1
VE1=scipy.io.loadmat('Results/QC/VE_PC1.mat')
VE1=VE1['VE1_PC1']
n_VE1_subs=VE1['AROMA1']; n_VE1_subs=n_VE1_subs[0,0]; n_VE1_subs=n_VE1_subs.shape[1]

curr_VE1=np.zeros([a_iters,n_VE1_subs])
mean_VE1=np.zeros([len(preprocs),n_VE1_subs])
min_VE1=np.zeros([len(preprocs),n_VE1_subs])
max_VE1=np.zeros([len(preprocs),n_VE1_subs])
for p in range(len(preprocs)):
        for a in range(10):
            curr=VE1['AROMA{}'.format(a+1)]
            curr=curr[0,0]
            curr_VE1[a,:]=curr[p,:]
        mean_vector = np.mean(curr_VE1, axis=1)
        mean_index = np.abs(mean_vector - np.mean(mean_vector)).argmin()
        mean_VE1[p,:]=curr_VE1[mean_index,:]
        min_index=np.argmin(mean_vector)
        min_VE1[p,:]=curr_VE1[min_index,:]
        max_index=np.argmax(mean_vector)
        max_VE1[p,:]=curr_VE1[max_index,:]
        
   

# Combine Results
mean_VE1=pd.DataFrame(mean_VE1)
mean_VE1.index=labels
mean_VE1=pd.DataFrame(mean_VE1.stack())
mean_VE1=mean_VE1.rename(columns={0: 'mean'})
mean_VE1['preprocs']=mean_VE1.index.get_level_values(0)

min_VE1=pd.DataFrame(min_VE1)
min_VE1=pd.DataFrame(min_VE1.stack())
min_VE1=min_VE1.rename(columns={0: 'min'})
max_VE1=pd.DataFrame(max_VE1)
max_VE1=pd.DataFrame(max_VE1.stack())
max_VE1=max_VE1.rename(columns={0: 'max'})

mean_VE1['min']=np.array(min_VE1['min']*100)
mean_VE1['max']=np.array(max_VE1['max']*100)
mean_VE1['mean']=mean_VE1['mean']*100

fig,ax=plt.subplots(1,1)
p=sns.boxenplot(data=mean_VE1, x='max', fill=False, y='preprocs', showfliers=False, palette=sns.color_palette(colors),ax=ax)
for artist in p.lines:
    artist.set_visible(False)
p=sns.boxenplot(data=mean_VE1, x='min', fill=False, y='preprocs', showfliers=False, palette=sns.color_palette(colors),ax=ax)
for artist in p.lines:
    artist.set_visible(False)
p2=sns.boxenplot(data=mean_VE1, x='mean', y='preprocs', showfliers=False, palette=sns.color_palette(colors), ax=ax)
#plt.title("Variance Explained by PC1")
plt.ylabel("")
plt.xlabel("Variance explained (%)")
plt.xlim([0,100])
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)



# Plot KRR across seeds
all_acc= scipy.io.loadmat('Results/KRR/age_accuracies.mat')
all_acc=all_acc['results']

dims=all_acc['AROMA1']
dims=dims[0,0]
dims=dims.shape
dim_vec=(dims[0]*dims[1] * dims[2])

curr_means=np.zeros([20,10])
mean_acc=np.zeros([len(preprocs),20])
min_acc=np.zeros([len(preprocs),20])
max_acc=np.zeros([len(preprocs),20])

for p in range(len(preprocs)):
    for a in range(10):
        curr=all_acc['AROMA{}'.format(a+1)]
        curr=curr[0,0]
        curr_proc=curr[:,:,:,p]
        curr_means_a=np.mean(curr_proc,axis=2)
        curr_means[:,a]=np.mean(curr_means_a,axis=1)
    
    mean_vector = np.mean(curr_means, axis=0)
    mean_index = np.abs(mean_vector - np.mean(mean_vector)).argmin()
    mean_acc[p,:]=curr_means[:,mean_index]
    min_index=np.argmin(mean_vector)
    min_acc[p,:]=curr_means[:,min_index]
    max_index=np.argmax(mean_vector)
    max_acc[p,:]=curr_means[:,max_index]

mean_acc=pd.DataFrame(mean_acc)
mean_acc.index=preprocs[:,0]
mean_acc=pd.DataFrame(mean_acc.stack())
mean_acc=mean_acc.rename(columns={0: 'mean'})
mean_acc['preprocs']=mean_acc.index.get_level_values(0)

min_acc=pd.DataFrame(min_acc)
min_acc=pd.DataFrame(min_acc.stack())
min_acc=min_acc.rename(columns={0: 'min'})
max_acc=pd.DataFrame(max_acc)
max_acc=pd.DataFrame(max_acc.stack())
max_acc=max_acc.rename(columns={0: 'max'})

mean_acc['min']=np.array(min_acc['min'])
mean_acc['max']=np.array(max_acc['max'])
mean_acc['mean']=mean_acc['mean']

fig,ax=plt.subplots(1,1)
p2=sns.boxplot(data=mean_acc, x='mean', y='preprocs', showfliers=False, palette=sns.color_palette(colors), ax=ax)

mean_means = mean_acc.groupby('preprocs')['mean'].mean()
min_means = mean_acc.groupby('preprocs')['min'].mean()
max_means = mean_acc.groupby('preprocs')['max'].mean()
categories = mean_means.index
# Manually draw means using scatter plot
ax.scatter(x=min_means, y=categories, color='darkblue', marker='^', s=20, label='Mean')
ax.scatter(x=max_means, y=categories, color='xkcd:gold', marker='^', s=20, label='Mean')
ax.scatter(x=mean_means, y=categories, color='white', marker='^', s=20, label='Mean')

plt.ylabel("")
plt.xlabel("mean accuracy (r)")
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.subplots_adjust(left=0.4)



##### Supplementary Figure 3
mean_mean=np.zeros(len(preprocs)-1)
max_mean=np.zeros(len(preprocs)-1)
min_mean=np.zeros(len(preprocs)-1)
for p in range(1,len(preprocs)):
    curr_name=str(preprocs[p])
    curr_name=curr_name.strip("[]").strip("'\"")
    
    curr_corrs = mmm_qc[mmm_qc['preprocs'] == curr_name]
    mean_mean[p-1]=np.mean(curr_corrs['mean'])
    max_mean[p-1]=np.mean(curr_corrs['max'])
    min_mean[p-1]=np.mean(curr_corrs['min'])
    

# Draw distributions of QC-FC.
# This is the old way - not using anymore.
# flat_corrs=mean_qc
# Li = -0.3  # lower index
# Ui = 0.3  # upper index
# r_sections = np.arange(Li, Ui + 0.025, 0.025)
# QC_distribution=np.zeros([len(preprocs), len(r_sections)])
# for k in range(len(preprocs)):
#     for r in range(len(r_sections) - 1):
#         vals_to_grab = np.logical_and(flat_corrs[k,:] >= r_sections[r], flat_corrs[k,:] <= r_sections[r + 1])
#         vals_in_range = flat_corrs[k,vals_to_grab]
#         QC_distribution[k, r] = np.sum(len(vals_in_range) / np.size(flat_corrs[k,:]) * 100)  # taking the % of all qcfc correlations

# ax=sns.heatmap(QC_distribution, cmap="Greens", xticklabels=np.round(r_sections,2), yticklabels=preprocs, 
#                cbar_kws={'label': '% of edges'}, linewidths=0.1, linecolor='grey')
# plt.title('HCP QC-FC Distributions')
# plt.xlabel('correlation strength (r)')


# # Plot KRR - also drawn an old wway
# all_acc= scipy.io.loadmat('Results/KRR/all_accuracies.mat')
# all_acc=all_acc['results']

# curr_means=np.zeros([a_iters,20])
# curr_m_iqr=np.zeros([a_iters,3])
# mean_acc=np.zeros([len(preprocs),3])
# min_acc=np.zeros([len(preprocs),3])
# max_acc=np.zeros([len(preprocs),3])
# for p in range(len(preprocs)):
#         for a in range(10):
#             curr=all_acc['AROMA{}'.format(a+1)]
#             curr=curr[0,0]
#             curr_proc=curr[:,:,p]
#             curr_means[a,:]=np.mean(curr_proc,axis=1)
#         curr_m_iqr[:,1],curr_m_iqr[:,2]=np.percentile(curr_means,[25,75],axis=1)
#         curr_m_iqr[:,0]=np.mean(curr_means,axis=1)
        
#         mean_acc[p,:]=np.mean(curr_m_iqr, axis=0)
#         min_acc[p,:]=np.min(curr_m_iqr, axis=0)
#         max_acc[p,:]=np.max(curr_m_iqr, axis=0)

# mean_acc=pd.DataFrame(mean_acc)
# mean_acc.index=preprocs[:,0]
# mean_acc=pd.DataFrame(mean_acc.stack())
# mean_acc=mean_acc.rename(columns={0: 'mean'})
# mean_acc['preprocs']=mean_acc.index.get_level_values(0)

# min_acc=pd.DataFrame(min_acc)
# min_acc=pd.DataFrame(min_acc.stack())
# min_acc=min_acc.rename(columns={0: 'min'})
# max_acc=pd.DataFrame(max_acc)
# max_acc=pd.DataFrame(max_acc.stack())
# max_acc=max_acc.rename(columns={0: 'max'})

# mean_acc['min']=np.array(min_acc['min'])
# mean_acc['max']=np.array(max_acc['max'])
# mean_acc['mean']=mean_acc['mean']

# fig,ax=plt.subplots(1,1)
# p=sns.boxenplot(data=mean_acc, x='preprocs', fill=False, y='max', showfliers=False, palette="tab10",ax=ax)
# p=sns.boxenplot(data=mean_acc, x='preprocs', fill=False, y='min', showfliers=False, palette="tab10",ax=ax)
# p2=sns.boxenplot(data=mean_acc, x='preprocs', y='mean', showfliers=False, palette="tab10", ax=ax)

# for artist in p.lines:
#     artist.set_visible(False)

# mean_means = mean_acc.groupby('preprocs')['mean'].mean()
# min_means = mean_acc.groupby('preprocs')['min'].mean()
# max_means = mean_acc.groupby('preprocs')['max'].mean()
# categories = mean_means.index
# # Manually draw means using scatter plot
# ax.scatter(x=categories, y=min_means, color='darkblue', marker='^', s=20, label='Mean')
# ax.scatter(x=categories, y=max_means, color='xkcd:gold', marker='^', s=20, label='Mean')
# ax.scatter(x=categories, y=mean_means, color='white', marker='^', s=20, label='Mean')

# plt.xlabel("")
# plt.xticks(rotation=70)
# plt.ylabel("mean accuracy (r)")
# plt.ylim([-0.04, 0.08])
# ax.spines['top'].set_visible(False)
# ax.spines['right'].set_visible(False)
