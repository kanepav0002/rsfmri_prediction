#!/bin/bash

# All the following is from the DiCER code, creates the restrictive GM map that DiCER uses to filter relevant voxels for clustering.
while getopts "a:f:l:s:" flag; do
	case "${flag}" in
		a) anat=${OPTARG} ;; 
		f) func=${OPTARG} ;;
		l) anat_loc=${OPTARG} ;;
		s) sub=${OPTARG} ;;
	esac
done


fast -o $anat_loc/$anat $anat_loc/$anat

base_anat="${anat%%.*}"

echo "Making Restrictive GM map for DiCER"
# Below is from DiCER code.
seg_temp=$anat_loc/"$base_anat"_seg.nii.gz
cp $anat_loc/"$base_anat"_pve_1.nii.gz $anat_loc/"$sub"_ds_gm.nii.gz
gm_mask_tmp=$anat_loc/"$sub"_gm_mask
fslmaths $anat_loc/"$sub"_ds_gm.nii.gz -thr 0.5 -bin $gm_mask_tmp

# Generate masks for GM to make it more restrictive:
mean_ts=$anat_loc/"$sub"_mean_ts
fslmaths $anat_loc/../func/$func -Tmean $mean_ts
# Taking the mean ts image, and just focusing in on grey matter
fslmaths $mean_ts -mul $gm_mask_tmp $mean_ts

# Now find the min/max
read min max <<< $(fslstats $mean_ts -r)

# Normalize the image and threshold the map to make a mask of epi of the top 60% of image intensity
gm_mask_restrictive=$anat_loc/"$sub"_mask_restrictive
fslmaths $mean_ts -div $max -thr 0.3 -bin $gm_mask_restrictive
fslmaths $gm_mask_restrictive -mul $gm_mask_tmp -mul 2 $gm_mask_restrictive

# Now we have everything to work with and combine it all together now (fix up too many volumes)
tissue_mask=$anat_loc/"$sub"_dtissue_func.nii.gz
fslmaths $seg_temp -add $gm_mask_restrictive $tissue_mask
fslroi $tissue_mask $tissue_mask 0 1

