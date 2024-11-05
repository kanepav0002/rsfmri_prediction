#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul  7 14:14:53 2023

@author: kanep
"""
import numpy as np
import pandas as pd
import os.path


sub_list=pd.read_csv('/home/kanep/kg98_scratch/Kane/UCLA/sub_list_full.csv', header=None)
sub_list=np.array(sub_list)

for idx in range(len(sub_list)):
    sub=str(sub_list[idx])
    end=len(sub)
    sub=sub[2:end-2]
  
    home_folder=('/home/kanep/kg98_scratch/Kane/UCLA/'+sub)
    
    if not os.path.isdir(home_folder):
        os.mkdir(home_folder)
        
    conf_folder=(home_folder+'/func/')
    confs=pd.read_csv((conf_folder+sub+'_task-rest_desc-confounds_timeseries.tsv'), sep='\t', engine='python')
    motion_params=confs[['rot_x', 'rot_y', 'rot_z','trans_x', 'trans_y', 'trans_z']]
    motion_params=motion_params.iloc[4:]
    fd=confs['framewise_displacement']
    fd=fd.iloc[4:]

    
    motion_params.to_csv(home_folder+'/func/prep_motion_params.par', header=None, index=False, sep='\t')
    fd.to_csv(home_folder+'/func/prep_fd.csv', header=None, index=False)
