#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  8 11:24:55 2023

@author: kanep
"""
from argparse import ArgumentParser
parser = ArgumentParser(epilog="resample_img.py -- Resample all imputs into same MNI space.")
parser.add_argument("-f", dest="func",
	help="fMRI file name", metavar="fMRI")
parser.add_argument("-r", dest="ref",
	help="Reference file", metavar="Standard Space")
parser.add_argument("-a", dest="anat",
	help="T1 Image", metavar="T1")

args = parser.parse_args()
# Setting the arguments
input_file = args.func
T1_file = args.anat
ref_file = args.ref

import nibabel as nib
from nilearn import image as nimg

func=nib.load(input_file)
T1=nib.load(T1_file)
ref=nib.load(ref_file)

resamp_func=nimg.resample_to_img(func,ref, interpolation='nearest')
resamp_t1=nimg.resample_to_img(T1, ref, interpolation='nearest')

out_T1=T1_file[:-7]
out_T1=out_T1+"_interp.nii.gz"
nib.save(resamp_t1,out_T1)

#out_func=input_file[:-7]
#out_func=out_func+"_interp.nii.gz"
#nib.save(resamp_func, out_func)

