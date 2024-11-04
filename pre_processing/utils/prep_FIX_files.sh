#!/bin/bash

while getopts "s:f:r:" flag; do
	case "${flag}" in
		s) sub=${OPTARG} ;; 
		f) func_loc=${OPTARG} ;; 
		r) repo_directory=${OPTARG} ;; 

	esac
done

mkdir -p $sub/reg
mkdir -p $sub/mc

cp $sub/func/ifr_func.nii.gz $sub/filtered_func_data.nii.gz

cd $sub
melodic -i filtered_func_data.nii.gz
cp filtered_func_data.ica/mean.nii.gz mean_func.nii.gz
cp filtered_func_data.ica/mask.nii.gz mask.nii.gz

# Get example function image (middle volume)
middle_vol=$(( $(fslinfo filtered_func_data.nii.gz | grep '^dim4' | awk '{print $2}') / 2 ))
fslroi filtered_func_data.nii.gz reg/example_func.nii.gz $middle_vol 1

# Copy structural image
cp anat/"$sub"_desc-preproc_T1w_brain.nii.gz reg/highres.nii.gz

# copy motion params
cp func/prep_motion_params.par mc/prefiltered_func_data_mcf.par

# copy in an identity matrix instead of affine transform, as the images are already in the same space
cp $repo_directory/dat/identity_matrix.mat reg/highres2example_func.mat
