#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  6 13:18:19 2024

@author: kanep
"""
import os
import numpy as np
import pandas as pd
import bct
from netneurotools import cluster
import scipy.io
from argparse import ArgumentParser
parser = ArgumentParser(epilog="mod_max_FC.py -- Does modularity maximisation on all FC .")
parser.add_argument("-f", dest="func_loc",
	help="FC location", metavar="func_loc")
parser.add_argument("-r", dest="reps",
	help="No. of DiCER repititions", metavar="reps")
parser.add_argument("-b", dest="base_name",
	help="FC file name", metavar="base_name")
parser.add_argument("-p", dest="parcels",
	help="No. of parcels", metavar="parcels")


args = parser.parse_args()
parcels=int(args.parcels)
reps=int(args.reps)
base_name=args.base_name

os.chdir(args.func_loc)

def consensus_Q(W,gamma,ci):
    W0 = W * (W > 0)
    s0 = np.sum(W0)
    B0 = W0 - gamma * np.outer(np.sum(W0, axis=1), np.sum(W0, axis=0)) / s0

    W1 = -W * (W < 0)
    s1 = np.sum(W1)
    if s1:
        B1 = W1 - gamma * np.outer(np.sum(W1, axis=1), np.sum(W1, axis=0)) / s1
    else:
        B1 = 0
    B = (B0 / s0) - (B1 / (s0 + s1))
    
    n = np.max(ci)
    b1 = np.zeros((n, n))
    for i in range(1, n + 1):
        for j in range(i, n + 1):
            # pool weights of nodes in same module
            bm = np.sum(B[np.ix_(ci == i, ci == j)])
            b1[i - 1, j - 1] = bm
            b1[j - 1, i - 1] = bm
    B = b1.copy()
    consensus_Q=np.trace(B)
    return consensus_Q

# Load initial data
All_FC=np.zeros([parcels,parcels,reps])
for e,i in enumerate(range(1,reps+1)):
    to_load=(base_name + str(i) + ".csv")
    tmp_FC=np.array(pd.read_csv(to_load, header=None))
    np.fill_diagonal(tmp_FC,1)
    All_FC[:,:,e]=tmp_FC
    
gamma=1.0
louv_reps=100
con_Q=np.zeros(reps)
All_Q=np.zeros([reps, louv_reps])
for p in range(reps):
    FC_s = All_FC[:, :, p]
    NGs_M=np.zeros([parcels,louv_reps])
    NGs_Q=np.zeros(louv_reps)
    for r in range(louv_reps):
        NGs_M[:, r], NGs_Q[r] = bct.algorithms.community_louvain(FC_s, gamma, B="negative_asym")
    ci = cluster.find_consensus(NGs_M)
    #num_communities_iter[p] = np.max(ci)
    con_Q[p] = consensus_Q(FC_s, gamma, ci)
    All_Q[p, :] = NGs_Q
    
indiv_max_Q, indiv_optimal_k = np.max(con_Q), np.argmax(con_Q)
#indiv_max_Q[1]=indiv_optimal_k
indiv_FC=All_FC[:,:,indiv_optimal_k]

indiv_optimal_k=indiv_optimal_k+1
np.savetxt('{}_opt.csv'.format(base_name), indiv_FC, delimiter=',')
print(indiv_optimal_k, file=open("opt_iters.txt", 'w'))
print(*con_Q, sep="," ,file=open("all_Q.txt", 'w'))

# Now do same for censored stuff
All_FC=np.zeros([parcels,parcels,reps])
for e,i in enumerate(range(1,reps+1)):
    to_load=(base_name + str(i) + "_cens.csv")
    tmp_FC=np.array(pd.read_csv(to_load, header=None))
    np.fill_diagonal(tmp_FC,1)
    All_FC[:,:,e]=tmp_FC
    
gamma=1.0
louv_reps=100
con_Q=np.zeros(reps)
All_Q=np.zeros([reps, louv_reps])
for p in range(reps):
    FC_s = All_FC[:, :, p]
    NGs_M=np.zeros([parcels,louv_reps])
    NGs_Q=np.zeros(louv_reps)
    for r in range(louv_reps):
        NGs_M[:, r], NGs_Q[r] = bct.algorithms.community_louvain(FC_s, gamma, B="negative_asym")
    ci = cluster.find_consensus(NGs_M)
    #num_communities_iter[p] = np.max(ci)
    con_Q[p] = consensus_Q(FC_s, gamma, ci)
    All_Q[p, :] = NGs_Q
    
indiv_max_Q, indiv_optimal_k = np.max(con_Q), np.argmax(con_Q)
#indiv_max_Q[1]=indiv_optimal_k
indiv_FC=All_FC[:,:,indiv_optimal_k]

indiv_optimal_k=indiv_optimal_k+1
np.savetxt('{}_opt_cens.csv'.format(base_name), indiv_FC, delimiter=',')
print(indiv_optimal_k, file=open("opt_iters_cens.txt", 'w'))
print(*con_Q, sep="," ,file=open("all_Q_cens.txt", 'w'))


