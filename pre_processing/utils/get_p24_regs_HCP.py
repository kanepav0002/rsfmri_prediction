#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 14 13:37:29 2024

@author: kanep
"""


import os
import pandas as pd
import numpy as np
from argparse import ArgumentParser
parser = ArgumentParser(epilog="This function collates some regressors for a single subject")
parser.add_argument("-s", dest="sub",
	help="Subject ID", metavar="sub")
parser.add_argument("-f", dest="folder",
	help="input folder", metavar="subject directory")
parser.add_argument("-m", dest="mov_params",
	help="movement parameters", metavar="movement parameters")
parser.add_argument("-w", dest="wm_file",
	help="white matter", metavar="white matter file")
parser.add_argument("-c", dest="csf_file",
	help="csf", metavar="csf file")
parser.add_argument("-g", dest="gs_file",
	help="global signal", metavar="global signal file")
parser.add_argument("-d", dest="dicer_file",
	help="dicer regressors", metavar="dicer regressors file")
args = parser.parse_args()
sub=args.sub
folder=args.folder

os.chdir(folder)
if os.path.exists(str(args.mov_params)):
    motion_params=np.array(pd.read_csv(args.mov_params, delimiter=r"\s+", header=None))
    motion_params=motion_params[:,0:6]
    p24_regs=motion_params
    for col in range(motion_params.shape[1]):
        curr_reg=motion_params[:,col]
        tmp_d1=np.zeros(curr_reg.shape[0])
        tmp_d1[1:]=np.diff(curr_reg)
    
        curr_sq=np.square(curr_reg)
        d1_sq=np.square(tmp_d1)
    
        p24_regs = np.column_stack((p24_regs, tmp_d1, curr_sq, d1_sq))
        out_regs=p24_regs
    
if os.path.exists(args.wm_file):
    WM=np.array(pd.read_csv(args.wm_file, header=None))
    WM_d1=np.zeros(WM.shape[0])
    WM_d1[1:]=np.diff(WM.T)

    if 'out_regs' in locals():
        out_regs=np.column_stack((out_regs, WM, WM_d1))
    else:
        out_regs=np.column_stack((WM, WM_d1))
        
if os.path.exists(args.csf_file):
    CSF=np.array(pd.read_csv(args.csf_file, header=None))
    CSF_d1=np.zeros(CSF.shape[0])
    CSF_d1[1:]=np.diff(CSF.T)
    
    if 'out_regs' in locals():
        out_regs=np.column_stack((out_regs, CSF, CSF_d1))
    else:
        out_regs=np.column_stack((CSF, CSF_d1))

if os.path.exists(args.gs_file):
    GS=np.array(pd.read_csv(args.gs_file, header=None))
    GS_d1=np.zeros(GS.shape[0])
    GS_d1[1:]=np.diff(GS.T)
    
    if 'out_regs' in locals():
        out_regs=np.column_stack((out_regs, GS, GS_d1))
    else:
        out_regs=np.column_stack((GS, GS_d1))
    
if os.path.exists(args.dicer_file):
    dic=np.array(pd.read_csv(args.dicer_file, delimiter=r"\s+", header=None))
    dic_d1=np.zeros(dic.shape)
    tmp=np.diff(dic.T)
    dic_d1[1:,:]=tmp.T
    
    if 'out_regs' in locals():
        out_regs=np.column_stack((out_regs, dic, dic_d1))
    else:
        out_regs=np.column_stack((dic, dic_d1))

np.savetxt('curr_regs.csv',out_regs,delimiter=',')

