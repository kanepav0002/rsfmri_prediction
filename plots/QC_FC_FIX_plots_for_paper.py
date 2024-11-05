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

os.chdir("/home/kanep/kg98_scratch/Kane/HCP")

############################################################################
# Legend
############################################################################
colors = [(1, 1, 0.3), (1, 0.95, 0.5), (0.176, 0.839, 0.341), (0.38, 0.969, 0.529), 
          (0.817, 0.188, 0.812), (0.961, 0.431, 0.914), (0.153, 0.729, 0.812), (0.396, 0.875, 0.941), 
          (0.902, 0.576, 0.094), (0.961, 0.773, 0.486), (0.878, 0.071, 0.145), (0.969, 0.38, 0.435), 
          (0.298, 0.31, 0.969), (0.392, 0.4, 0.961),(0.678, 0.659, 0.678), (0.839, 0.824, 0.839),
          (0.96, 0.15, 0.52), (1, 0.41, 0.7)]
labels = [
    "1. 24P", "2. 28P", "3. FIX/AROMA + 28P", "4. FIX/AROMA + 28P censored",
    "5. FIX/AROMA + 28P + GSR", "6. FIX/AROMA + 28P + GSR censored",
    "7. FIX/AROMA + 28P + 1DiCER", "8. FIX/AROMA + 28P + 1DiCER censored",
    "9. FIX/AROMA + 28P + 2DiCER", "10. FIX/AROMA + 28P + 2DiCER censored",
    "11. FIX/AROMA + 28P + 3DiCER", "12. FIX/AROMA + 28P + 3DiCER censored",
    "13. FIX/AROMA + 28P + 4DiCER", "14. FIX/AROMA + 28P + 4DiCER censored",
    "15. FIX/AROMA + 28P + 5DiCER", "16. FIX/AROMA + 28P + 5DiCER censored",
    "17. FIX/AROMA + 28P + t-DiCER", "18. FIX/AROMA + 28P + t-DiCER censored"
]

# Create a figure
plt.figure(figsize=(6, 8))

# Create the legend manually
for i, label in enumerate(labels):
    plt.scatter([], [], color=colors[i], label=label, s=100, marker='s')

plt.legend(frameon=False, loc='center', fontsize=10)
plt.axis('off')
plt.show()

# Figure 1
GSP_FD=pd.read_csv('/home/kanep/kg98_scratch/Kane/GSP/GSP_meanFD.csv', header=None)
GSP_FD=GSP_FD.T
HCP_FD=pd.read_csv('/home/kanep/kg98_scratch/Kane/HCP/HCP_meanFD.csv', header=None)
HCP_FD=HCP_FD.T
CNP_FD=pd.read_csv('/home/kanep/kg98_scratch/Kane/UCLA/CNP_meanFD.csv', header=None)
CNP_FD=CNP_FD.T

gsp_data = GSP_FD[0]
hcp_data = HCP_FD[0]
cnp_data = CNP_FD[0]

fig, axs = plt.subplots(1, 3, figsize=(18, 5), sharey=True)
# Plot GSP data
sns.kdeplot(gsp_data, ax=axs[0], shade=True, color='blue')
axs[0].set_title('GSP')
axs[0].set_ylabel('Kernel Density Estimate')
axs[0].set_xlabel('')
axs[0].axvline(x=0.3, color='red', linestyle='--')
axs[0].legend()

# Plot HCP data
sns.kdeplot(hcp_data, ax=axs[1], shade=True, color='purple')
axs[1].set_title('HCP')
axs[1].set_xlabel('Mean FD')
axs[1].axvline(x=0.3, color='red', linestyle='--')
axs[1].legend()

# Plot CNP data
sns.kdeplot(cnp_data, ax=axs[2], shade=True, color='green')
axs[2].set_title('CNP')
axs[2].set_xlabel('')
axs[2].axvline(x=0.3, color='red', linestyle='--', label='cutoff')
axs[2].legend()

# Remove top and right borders
for ax in axs:
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

plt.tight_layout()
plt.show()


###########################
# QC - FC
##############################
colors = [(1, 1, 0.3), (1, 0.95, 0.5), (0.176, 0.839, 0.341), (0.38, 0.969, 0.529), 
          (0.817, 0.188, 0.812), (0.961, 0.431, 0.914), (0.153, 0.729, 0.812), (0.396, 0.875, 0.941), 
          (0.902, 0.576, 0.094), (0.961, 0.773, 0.486), (0.878, 0.071, 0.145), (0.969, 0.38, 0.435), 
          (0.298, 0.31, 0.969), (0.392, 0.4, 0.961),(0.678, 0.659, 0.678), (0.839, 0.824, 0.839),
          (0.96, 0.15, 0.52), (1, 0.41, 0.7)]

labels = np.array([("1"), ("2"), ("3**"), ("4**"), ("5**"), ("6**"), ("7**"), ("8**"), ("9**"), ("10**"),
          ("11**"), ("12**"), ("13**"), ("14**"), ("15**"), ("16**"), ("17**"), ("18**")])
# Load results
QC_corrs=scipy.io.loadmat('Results/QC/QC_FC_corrs.mat')
QC_corrs=QC_corrs['QC_FC']
preprocs=np.array(pd.read_csv("Results/QC/preprocs.csv", header=None))
n_edges=QC_corrs; n_edges=n_edges.shape[1]

mean_qc=QC_corrs
mean_qc=pd.DataFrame(mean_qc)
mean_qc.index=labels #preprocs[:,0] # Change this to re-label.
mean_qc=pd.DataFrame(mean_qc.stack())
mean_qc=mean_qc.rename(columns={0: 'mean'})
mean_qc['preprocs']=mean_qc.index.get_level_values(0)

# Plot the correlation distributions
sns.set_theme(style="white", rc={"axes.facecolor": (0, 0, 0, 0), 'axes.linewidth':2})
palette=sns.color_palette(colors)
g = sns.FacetGrid(mean_qc, palette=palette, row='preprocs', hue='preprocs', aspect=6, height=2.2)
g.map_dataframe(sns.kdeplot, x='mean', fill=True, alpha=1)
g.map_dataframe(sns.kdeplot, x='mean', color='black')
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
plt.suptitle("", fontsize=20, x=0.2, y=0.95)
plt.subplots_adjust(bottom=0.15)

# Load Distance Dependence
# This block has to be run with seaborn version 0.11 (module load anaconda 2019 version)
# The xerr doesn't play nicely with later versions.
dist_dep=scipy.io.loadmat('Results/QC/distance_dependence.mat')
dist_dep=dist_dep['dist_dep']

mmm_distdep=pd.DataFrame()
mmm_distdep['preprocs']=labels
mmm_distdep['mean']=dist_dep[0]

plt.figure()
ax=sns.barplot(data=mmm_distdep, x='mean', y='preprocs', palette=sns.color_palette(colors))
#ax.bar_label(ax.containers[0], fmt='%.3f')
plt.xlim(-0.3,0.2)
plt.axvline(x=0, color='grey', linestyle='-')
plt.xlabel('Spearmans rho (r)')
plt.ylabel=(" ")
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.subplots_adjust(left=0.4)

# Draw VE by PC1
VE1=pd.read_csv('Results/VE1/VE1_PC1.csv',header=None)
#VE1=VE1[1:]
VE1.index=labels
VE1_df= pd.DataFrame(VE1.stack())

VE1_df = VE1_df.rename(columns={0: 'Variance'})
VE1_df["Pipeline_v"] = VE1_df.index.get_level_values(0)
VE1_df["Variance"]=VE1_df["Variance"]*100

sns.catplot(data=VE1_df, x='Variance', y='Pipeline_v', kind='boxen', palette=sns.color_palette(colors), legend=None, showfliers=False)
#plt.ylabel("")
plt.xlabel("Variance explained (%)")
plt.xlim([0,100])


# Plot KRR
# all_acc= scipy.io.loadmat('Results/KRR/all_accuracies.mat')
# all_acc=all_acc['results']

# curr_m_iqr=np.zeros([1,3])
# mean_acc=np.zeros([len(preprocs),3])
# for p in range(len(preprocs)):
#     #curr=all_acc['FIX']
#     #curr=curr[0,0]
#     curr_proc=all_acc[:,:,p]
#     curr_means=np.mean(curr_proc,axis=1)
#     curr_m_iqr[:,1],curr_m_iqr[:,2]=np.percentile(curr_means,[25,75],axis=0)
#     curr_m_iqr[:,0]=np.mean(curr_means,axis=0)
        
#     mean_acc[p,:]=np.mean(curr_m_iqr, axis=0)


# mean_acc=pd.DataFrame(mean_acc)
# mean_acc.index=preprocs[:,0]
# mean_acc=pd.DataFrame(mean_acc.stack())
# mean_acc=mean_acc.rename(columns={0: 'mean'})
# mean_acc['preprocs']=mean_acc.index.get_level_values(0)
# mean_acc['mean']=mean_acc['mean']

# fig,ax=plt.subplots(1,1)
# p2=sns.boxplot(data=mean_acc, x='preprocs', y='mean', showfliers=False, palette="tab10", ax=ax)
# #for artist in p2.lines:
# #    artist.set_visible(False)
# mean_means = mean_acc.groupby('preprocs')['mean'].mean()
# categories = mean_means.index
# plt.scatter(x=categories, y=mean_means, color='green', marker='^', s=20, label='Mean')


# plt.xlabel("")
# plt.xticks(rotation=90)
# plt.ylabel("mean accuracy (r)")
# ax.spines['top'].set_visible(False)
# ax.spines['right'].set_visible(False)

# Plot KRR across seeds
all_acc= scipy.io.loadmat('Results/KRR/all_accuracies.mat')
all_acc=all_acc['results']

cog_index=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,24,25,29]
personality_index=[30,31,32,33,34,41,42,43,44,45,46,49]

dims=all_acc.shape
dim_vec=(dims[0]*dims[1] * dims[2])
to_plot=np.zeros([20, len(preprocs)])
for p in range(len(preprocs)):
    curr_proc=all_acc[:,:,:,p]
    curr_proc=curr_proc[:,cog_index,:] # Remove if doing all behaviours
    curr_means=np.mean(curr_proc,axis=2)
    curr_means=np.mean(curr_means,axis=1)
    
    to_plot[:,p] = curr_means
    
mean_acc=pd.DataFrame(to_plot)
mean_acc.columns=labels
mean_acc=mean_acc.stack().reset_index(drop=False)
mean_acc=mean_acc.drop(columns="level_0")
mean_acc.columns=["preprocs", "mean"]


fig,ax=plt.subplots(1,1)
p2=sns.boxplot(data=mean_acc, x='mean', y='preprocs', showfliers=False, palette=sns.color_palette(colors), ax=ax)

#for artist in p2.lines:
#    artist.set_visible(False)
mean_means = mean_acc.groupby('preprocs')['mean'].mean()
categories = mean_means.index
plt.scatter(x=mean_means, y=categories, color='white', marker='^', s=20, label='Mean', zorder=5)


#plt.ylabel("")
plt.xlabel("mean accuracy (r)")
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.subplots_adjust(left=0.4)
