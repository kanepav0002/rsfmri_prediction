#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 25 10:18:17 2024

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
parser.add_argument("-d", dest="dicer_file",
	help="dicer regressors", metavar="dicer regressors file")
args = parser.parse_args()
sub=args.sub
folder=args.folder

os.chdir(folder)

dic=np.array(pd.read_csv(args.dicer_file, delimiter=r"\s+", header=None))
#dic_d1=np.zeros(dic.shape)
#tmp=np.diff(dic.T)
#dic_d1[1:,:]=tmp.T
    
np.savetxt('Dic_regs.csv',dic,delimiter=',')

