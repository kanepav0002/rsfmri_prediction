#!/bin/bash

while getopts "a:f:r:" flag; do
	case "${flag}" in
		a) T1=${OPTARG} ;;
		f) func_loc=${OPTARG} ;;
		r) repo_directory=${OPTARG} ;;
	esac
done


anat_loc=${func_loc::-4}
anat_loc=$anat_loc"anat"

matlab -nodisplay -r "segment_anat('$anat_loc/$T1".gz"', '$repo_directory'); exit;"
cd $anat_loc
gm="c1"$T1
wm="c2"$T1
csf="c3"$T1

# Here we are getting ventricle CSF as the regressor
#Get binary mask from brain extracted T1 and from thresholded gm
fslmaths $T1 -bin brain_mask
fslmaths $gm -thr 0.95 -bin vmask

# Dilate gm mask twice
fslmaths vmask -dilD vmask
fslmaths vmask -dilD vmask
# Add WM and invert
fslmaths vmask -add $wm -binv vmask

# Erode whole brain mask twice
fslmaths brain_mask -eroF e_brain_mask
fslmaths brain_mask -eroF e_brain_mask

csf="v"$csf
# multiply eroded brain mask with vmask
fslmaths vmask -mul e_brain_mask $csf

# erode csf twice
fslmaths $csf -eroF -bin csf_e1
fslmaths csf_e1 -eroF csf_e2

# erode wm 5 times
fslmaths $wm -eroF -bin wm_e1
fslmaths wm_e1 -eroF -bin wm_e2
fslmaths wm_e2 -eroF -bin wm_e3
fslmaths wm_e3 -eroF -bin wm_e4
fslmaths wm_e4 -eroF -bin wm_e5

# check there is at least 5 voxels left in the csf / wm masks or take and ealier erosion iteration.
cd $repo_directory/pre_processing/utils
csf=$csf".gz"
matlab -nodisplay -r "check_erosions('$anat_loc', '$csf'); exit;"

# Delete intermediate files
cd $anat_loc
rm wm_e?.nii*
rm csf_e?.nii*







