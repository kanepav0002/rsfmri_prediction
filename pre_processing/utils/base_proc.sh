#!/bin/bash

while getopts "s:f:a:p:r:n:" flag; do
	case "${flag}" in
		s) sub=${OPTARG} ;; 
		f) func_loc=${OPTARG} ;; 
		a) anat_loc=${OPTARG} ;; 
		p) parc_file=${OPTARG} ;; 
		r) repo_directory=${OPTARG} ;;
		n) num_frames=${OPTARG} ;;
	esac
done

cd $repo_directory/pre_processing/utils

# Remove first 4 frames
to_rmv=$((num_frames-4))
fslroi $func_loc/"$sub"_task-rest_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz $func_loc/fr_func.nii.gz 4 $to_rmv

# remove skull from anat using fmriprep masks - often cuts too much brain off prefer to use bet
#fslmaths "$anat_loc/"$sub"_desc-preproc_T1w.nii.gz" -mas $anat_loc/sub-10159_desc-brain_mask.nii.gz "$anat_loc/"$sub"_desc-preproc_T1w_brain.nii.gz"

# Remove skull from anat with bet
bet $anat_loc/"$sub"_desc-preproc_T1w.nii.gz $anat_loc/"$sub"_desc-preproc_T1w_brain.nii.gz

T1="$anat_loc/"$sub"_desc-preproc_T1w_brain.nii.gz"
func_file="fr_func"
rm $func_loc/fr_func.nii

# Put anat in MNI func space 
flirt -in $func_loc/$func_file".nii.gz" -ref $T1  -omat $anat_loc/"$sub"_func2struct.mat -dof 12
flirt -ref ${FSLDIR}/data/standard/MNI152_T1_2mm_brain -in $T1 -omat $anat_loc/"$sub"_affine_tranf.mat -dof 12

fnirt --in=$anat_loc/"$sub"_desc-preproc_T1w_brain.nii.gz --aff=$anat_loc/"$sub"_affine_tranf.mat --cout=$anat_loc/"$sub"_nonlinear_transf --config=T1_2_MNI152_2mm --iout=$anat_loc/"$sub"_desc-preproc_T1w_brain.nii.gz
applywarp --ref=${FSLDIR}/data/standard/MNI152_T1_2mm_brain --in=$func_loc/$func_file".nii.gz" --warp=$anat_loc/"$sub"_nonlinear_transf --premat=$anat_loc/"$sub"_func2struct.mat --out=$func_loc/$func_file".nii.gz"

# remove skull from func
fslmaths $T1 -bin $func_loc/../anat/bmask.nii.gz
fslmaths $func_loc/$func_file -mas $func_loc/../anat/bmask.nii.gz $func_loc/$func_file

# Intensity Normalise
matlab -nodisplay -r "Inorm('$func_file".nii.gz"', '$func_loc', '$repo_directory'); exit;"

# Get WM / CSF / GM masks
T1_name="$sub"_desc-preproc_T1w_brain.nii
sh linden_WM_CSF_masks.sh -a $T1_name -f $func_loc -r $repo_directory
rm $func_loc/ifr_func.nii $func_loc/ifr_func.mat
# Create grey matter probability for parcellation
gm="$anat_loc/c1"$sub"_desc-preproc_T1w_brain"
fslmeants -i "$gm".nii.gz --label=$parc_file -o $func_loc/fsl_gm_prob.txt

# Get Restrictive grey matter mask for DiCER
sh Make_Restrictive_GM_Map.sh -a "$T1_name".gz -f ifr_func.nii.gz -l $anat_loc -s $sub
