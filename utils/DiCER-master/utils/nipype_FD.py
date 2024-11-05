#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul  7 14:14:53 2023

@author: kanep
"""

from nipype.algorithms.confounds import FramewiseDisplacement
import pandas as pd

sub='sub-10159'
home_folder=('/home/kanep/kg98_scratch/Kane/GSP/'+sub+'/func/')
func_file=('/home/kanep/kg98_scratch/kevo/UCLA_downloaded/ds000030_R1.0.5/'+sub+'/func/'+sub+'_task-rest_bold.nii.gz')
confs=pd.read_csv((home_folder+sub+'_task-rest_desc-confounds_timeseries.tsv'), sep='\t', engine='python')
motion_params=confs[['rot_x', 'rot_y', 'rot_z','trans_x', 'trans_y', 'trans_z']]
motion_params.to_csv(home_folder+'prep_motion_params.txt', header=None, index=False, sep='\t')


sub_list=pd.read_csv('/home/kanep/kg98_scratch/Kane/UCLA/sub_list_full.csv', header=None)
for s in range(len(sub_list)):
    sub=sub_list.iloc[s,0]
    home_folder=('/home/kanep/kg98_scratch/Kane/UCLA/'+sub+'/func/')
    confs=pd.read_csv((home_folder+sub+'_task-rest_desc-confounds_timeseries.tsv'), sep='\t', engine='python')
    motion_params=confs[['rot_x', 'rot_y', 'rot_z','trans_x', 'trans_y', 'trans_z']]
    motion_params.to_csv(home_folder+'prep_motion_params.par', header=None, index=False, sep='\t')

    
    # home_folder=('/home/kanep/kg98_scratch/Kane/UCLA/'+sub+'/func/')
    # fd=FramewiseDisplacement(in_file=('/home/kanep/kg98/Kane/DiCER/UCLA_proc/'+sub+'/orig_mcf.par'), parameter_source='FSL',out_file=(home_folder+'nipype_mcflirt_fd.txt'))
    # fd.run()
    
    # nip_fd=pd.read_csv(home_folder+'nipype_mcflirt_fd.txt', header=None)
    # nip_fd.iloc[0,0]=0
    # nip_fd.to_csv(home_folder+'nipype_mcflirt_fd.txt', header=None, index=False)
print('done')