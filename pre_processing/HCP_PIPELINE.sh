#!/bin/bash

#SBATCH --job-name=HCP_PIPELINE
#SBATCH --time=0-40:30:00
#SBATCH --ntasks=1
#SBATCH --mem=75000
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

while getopts "s:r:" flag; do
	case "${flag}" in
		s) sub=${OPTARG} ;; 
		r) reps=${OPTARG} ;;
	esac
done

# Manual inputs
repo_directory='location the repo has been downloaded to'
func_loc="subjects functional directory e.g. /home/HCP/$sub/func"
anat_loc="subjects anatomical directory e.g. /home/HCP/$sub/anat"
TR="TR of the fmri scan e.g. 3"
TR_ms="TR of the fmri scan in ms"

# Modules required - make sure these exist
# Replace these with the module names in your cluster
module purge 
module load fsl/6.0.1
module load matlab

#activate python
source $repo_directory/DiCER_env/bin/activate

dicer_folder=$repo_directory/utils/DiCER-master/

cd $repo_directory/pre_processing/utils

parc_file="$repo_directory/dat/Schaefer2018_300Parcels_17Networks_order_FSLMNI152_2mm.nii.gz"
num_parcels=300
func_file="rfMRI_REST1_LR.nii.gz"


# Downsample T1 from 1mm to 2mm
python resample_images.py -f $func_loc/rfMRI_REST1_LR.nii.gz -r $repo_directory/dat/MNI152_T1_2mm_brain.nii.gz -a $anat_loc/T1w_brain.nii.gz

# Make brain mask & WM/CSF/GM masks & GM map
fslmaths $anat_loc/T1w_brain_interp.nii.gz -bin $anat_loc/brain_mask.nii.gz
brain_mask=$anat_loc/brain_mask.nii.gz
T1_name="T1w_brain_interp.nii"
sh linden_WM_CSF_masks.sh -a $T1_name -f $func_loc -r $repo_directory
# Create grey matter probability for parcellation
gm="$anat_loc/c1T1w_brain_interp"
fslmeants -i "$gm".nii.gz --label=$parc_file -o $func_loc/fsl_gm_prob.txt

##############################################################
# Base Processing Pre-FIX
# Will do 24P regression w/ and w/out 2P regression and save FC
#############################################################

fslmeants -i $func_loc/rfMRI_REST1_LR.nii.gz -o $func_loc/WM_pre_fix_ts.txt -m $anat_loc/final_wmmask.nii.gz
fslmeants -i $func_loc/rfMRI_REST1_LR.nii.gz -o $func_loc/CSF_pre_fix_ts.txt -m $anat_loc/final_csfmask.nii.gz

python get_p24_regs_HCP.py -s $sub -f $func_loc -m $func_loc/Movement_Regressors.txt -w $func_loc/WM_pre_fix_ts.txt -c $func_loc/CSF_pre_fix_ts.txt -g "not_this_time" -d "not_this_time"

# Create a list with no censored frames
dim_time=$(fslinfo $func_loc/rfMRI_REST1_LR.nii.gz | awk '/dim4/ {print $2}')
dim_time=${dim_time:0:4}
rm $func_loc/no_censor.txt
for ((i=0; i<dim_time; i++)); do
    echo "1" >> $func_loc/no_censor.txt
done


# Do Regression for 24p 
python get_p24_regs_HCP.py -s $sub -f $func_loc -m $func_loc/Movement_Regressors.txt -w "not_this_time" -c "not_this_time" -g "not_this_time" -d "not_this_time"
echo $func_loc/rfMRI_REST1_LR.nii.gz > $func_loc/tmp_fmri_list.txt
echo $func_loc/24P.nii.gz > $func_loc/tmp_out_list.txt
echo $func_loc/curr_regs.csv > $func_loc/tmp_regs.txt
echo $func_loc/no_censor.txt > $func_loc/tmp_c_frames_list.txt
matlab -nodisplay -r "CBIG_glm_regress_vol('$func_loc/tmp_fmri_list.txt', '$func_loc/tmp_out_list.txt','$func_loc/tmp_regs.txt', '1', '$func_loc/tmp_c_frames_list.txt', '1', '0', '$repo_directory'); exit;"
# And 24 + 2P
python get_p24_regs_HCP.py -s $sub -f $func_loc -m $func_loc/Movement_Regressors.txt -w $func_loc/WM_pre_fix_ts.txt -c $func_loc/CSF_pre_fix_ts.txt -g "not_this_time" -d "not_this_time"
echo $func_loc/rfMRI_REST1_LR.nii.gz > $func_loc/tmp_fmri_list.txt
echo $func_loc/24_2P.nii.gz > $func_loc/tmp_out_list.txt
echo $func_loc/curr_regs.csv > $func_loc/tmp_regs.txt
echo $func_loc/no_censor.txt > $func_loc/tmp_c_frames_list.txt
matlab -nodisplay -r "CBIG_glm_regress_vol('$func_loc/tmp_fmri_list.txt', '$func_loc/tmp_out_list.txt','$func_loc/tmp_regs.txt', '1', '$func_loc/tmp_c_frames_list.txt', '1', '0', '$repo_directory'); exit;"

# Bandpass Filter 
matlab -nodisplay -r "bandpass_filter('$func_loc', '24P.nii.gz', '$brain_mask', '$TR', '0.008', '0.08', '$repo_directory'); exit;"
matlab -nodisplay -r "bandpass_filter('$func_loc', '24_2P.nii.gz', '$brain_mask', '$TR', '0.008', '0.08', '$repo_directory'); exit;"

# Make FC from those two files
fslmeants -i "$gm".nii.gz --label=$parc_file -o $func_loc/fsl_gm_prob.txt
fslmaths $func_loc/24P_bpass.nii.gz -mul "$gm".nii.gz $func_loc/epi_w
fslmeants -i $func_loc/epi_w --label=$parc_file -o $func_loc/fsl_parc_ts.txt
matlab -nodisplay -r "linden_parc('$func_loc', 'Schaef_300_24P.csv'); exit;"

fslmaths $func_loc/24_2P_bpass.nii.gz -mul "$gm".nii.gz $func_loc/epi_w
fslmeants -i $func_loc/epi_w --label=$parc_file -o $func_loc/fsl_parc_ts.txt
matlab -nodisplay -r "linden_parc('$func_loc', 'Schaef_300_24P_2P.csv'); exit;"


##############################################################
# Do FIX on 24P data
# do one pipeline with Just FIX
# Another where an additional regression of 24p and 2p.
# Split into scrub and no scrub.
# Another with 24p, 2p, GSR
# Another with 24p 2p, 1DiCER
#############################################################


# Get new regs
fslmeants -i $func_loc/rfMRI_REST1_LR_hp2000_clean.nii.gz -o $func_loc/WM_post_fix_ts.txt -m $anat_loc/final_wmmask.nii.gz
fslmeants -i $func_loc/rfMRI_REST1_LR_hp2000_clean.nii.gz -o $func_loc/CSF_post_fix_ts.txt -m $anat_loc/final_csfmask.nii.gz
fslmeants -i $func_loc/rfMRI_REST1_LR_hp2000_clean.nii.gz -o $func_loc/GS_post_fix_ts.txt -m $anat_loc/brain_mask.nii.gz


# Do Regressions with and without ignoring scrubbed frames
python get_p24_regs_HCP.py -s $sub -f $func_loc -m "not_this_time" -w $func_loc/WM_post_fix_ts.txt -c $func_loc/CSF_post_fix_ts.txt -g "not_this_time" -d "not_this_time"
echo $func_loc/rfMRI_REST1_LR_hp2000_clean.nii.gz > $func_loc/tmp_fmri_list.txt
echo $func_loc/FIX_24P_2P.nii.gz > $func_loc/tmp_out_list.txt
echo $func_loc/curr_regs.csv > $func_loc/tmp_regs.txt
echo $func_loc/no_censor.txt > $func_loc/tmp_c_frames_list.txt
matlab -nodisplay -r "CBIG_glm_regress_vol('$func_loc/tmp_fmri_list.txt', '$func_loc/tmp_out_list.txt','$func_loc/tmp_regs.txt', '1', '$func_loc/tmp_c_frames_list.txt', '1', '0', '$repo_directory'); exit;"
echo $func_loc/FIX_24P_2P_c.nii.gz > $func_loc/tmp_out_list.txt
echo $func_loc/censored_frames.txt > $func_loc/tmp_c_frames_list.txt
matlab -nodisplay -r "CBIG_glm_regress_vol('$func_loc/tmp_fmri_list.txt', '$func_loc/tmp_out_list.txt','$func_loc/tmp_regs.txt', '1', '$func_loc/tmp_c_frames_list.txt', '1', '0', '$repo_directory'); exit;"

python get_p24_regs_HCP.py -s $sub -f $func_loc -m "not_this_time" -w $func_loc/WM_post_fix_ts.txt -c $func_loc/CSF_post_fix_ts.txt -g $func_loc/GS_post_fix_ts.txt -d "not_this_time"
echo $func_loc/rfMRI_REST1_LR_hp2000_clean.nii.gz > $func_loc/tmp_fmri_list.txt
echo $func_loc/FIX_24P_2P_GS.nii.gz > $func_loc/tmp_out_list.txt
echo $func_loc/curr_regs.csv > $func_loc/tmp_regs.txt
echo $func_loc/no_censor.txt > $func_loc/tmp_c_frames_list.txt
matlab -nodisplay -r "CBIG_glm_regress_vol('$func_loc/tmp_fmri_list.txt', '$func_loc/tmp_out_list.txt','$func_loc/tmp_regs.txt', '1', '$func_loc/tmp_c_frames_list.txt', '1', '0', '$repo_directory'); exit;"
echo $func_loc/FIX_24P_2P_GS_c.nii.gz > $func_loc/tmp_out_list.txt
echo $func_loc/censored_frames.txt > $func_loc/tmp_c_frames_list.txt
matlab -nodisplay -r "CBIG_glm_regress_vol('$func_loc/tmp_fmri_list.txt', '$func_loc/tmp_out_list.txt','$func_loc/tmp_regs.txt', '1', '$func_loc/tmp_c_frames_list.txt', '1', '0', '$repo_directory'); exit;"

# Scrub the data
matlab -nodisplay -r "CBIG_preproc_censor_wrapper('$func_loc/FIX_24P_2P_GS_c.nii.gz', '$func_loc/censored_frames.txt', '$TR_ms', '$func_loc/tmp_interp_bold.nii.gz', '$func_loc/FIX_24P_2p_GS_cens.nii.gz', '$brain_mask', '10', '$repo_directory'); exit;"
matlab -nodisplay -r "CBIG_preproc_censor_wrapper('$func_loc/FIX_24P_2P_c.nii.gz', '$func_loc/censored_frames.txt', '$TR_ms', '$func_loc/tmp_interp_bold.nii.gz', '$func_loc/FIX_24P_2P_cens.nii.gz', '$brain_mask', '10', '$repo_directory'); exit;"

# Do DiCER
cd $dicer_folder
echo "Starting DiCER"
sh DiCER_very_lightweight.sh -i FIX_24P_2P.nii.gz  -a ../anat/T1w_brain_interp.nii.gz -w $func_loc -s $sub -r $reps -e 1.0

# Move the d_tissue file to anat 
mv $func_loc/"$sub"_dtissue_func.nii.gz $anat_loc/"$sub"_dtissue_func.nii.gz

# Collate Regressors and get first derivatives.
cd $repo_directory/pre_processing/utils

for dic_num in $(seq 1 $reps); do

	python get_dic_regs.py -s $sub -f $func_loc -d $func_loc/"$sub"_dbscan_"$dic_num"regressors.tsv

	echo $func_loc/FIX_24P_2P.nii.gz > $func_loc/tmp_fmri_list.txt
	echo $func_loc/FIX_24P_2P_Dic"$dic_num".nii.gz > $func_loc/tmp_out_list.txt
	echo $func_loc/Dic_regs.csv > $func_loc/tmp_regs.txt
	matlab -nodisplay -r "CBIG_glm_regress_vol('$func_loc/tmp_fmri_list.txt', '$func_loc/tmp_out_list.txt','$func_loc/tmp_regs.txt', '1', '$func_loc/tmp_c_frames_list.txt', '1', '0', '$repo_directory'); exit;"
	echo $func_loc/FIX_24P_2P_Dic"$dic_num"_c.nii.gz > $func_loc/tmp_out_list.txt
	echo $func_loc/censored_frames.txt > $func_loc/tmp_c_frames_list.txt
	matlab -nodisplay -r "CBIG_glm_regress_vol('$func_loc/tmp_fmri_list.txt', '$func_loc/tmp_out_list.txt','$func_loc/tmp_regs.txt', '1', '$func_loc/tmp_c_frames_list.txt', '1', '0', '$repo_directory'); exit;"

	matlab -nodisplay -r "CBIG_preproc_censor_wrapper('$func_loc/FIX_24P_2P_Dic"$dic_num"_c.nii.gz', '$func_loc/censored_frames.txt', '$TR_ms', '$func_loc/tmp_interp_bold.nii.gz', '$func_loc/FIX_24P_2P_Dic"$dic_num"_cens.nii.gz', '$brain_mask', '10', '$repo_directory'); exit;"

	matlab -nodisplay -r "bandpass_filter('$func_loc', 'FIX_24P_2P_Dic"$dic_num".nii.gz', '$brain_mask', '$TR', '0.008', '0.08', '$repo_directory'); exit;"
	matlab -nodisplay -r "bandpass_filter('$func_loc', 'FIX_24P_2P_Dic"$dic_num"_cens.nii.gz', '$brain_mask', '$TR', '0.008', '0.08', '$repo_directory'); exit;"
	
	fslmaths $func_loc/FIX_24P_2P_Dic"$dic_num"_bpass.nii.gz -mul "$gm".nii.gz $func_loc/epi_w
	fslmeants -i $func_loc/epi_w --label=$parc_file -o $func_loc/fsl_parc_ts.txt
	matlab -nodisplay -r "linden_parc('$func_loc', 'Schaef_300_FIX_24P_2P_Dic"$dic_num".csv'); exit;"

	fslmaths $func_loc/FIX_24P_2P_Dic"$dic_num"_cens_bpass.nii.gz -mul "$gm".nii.gz $func_loc/epi_w
	fslmeants -i $func_loc/epi_w --label=$parc_file -o $func_loc/fsl_parc_ts.txt
	matlab -nodisplay -r "linden_parc('$func_loc', 'Schaef_300_FIX_24P_2P_Dic"$dic_num"_cens.csv'); exit;"
done

# Get modularity Maximised FC
python mod_max_FC.py -f $func_loc -r $reps -b "FC_Schaef_300_FIX_24P_2P_Dic" -p 300

# bandpass filter and parcellate all FIX files
fmri_files=("rfMRI_REST1_LR_hp2000_clean" "FIX_24P_2P" "FIX_24P_2P_cens" "FIX_24P_2P_GS" "FIX_24P_2p_GS_cens")
FC_files=("Schaef_300_FIX.csv" "Schaef_300_FIX_24P_2P.csv" "Schaef_300_FIX_24P_2P_cens.csv" "Schaef_300_FIX_24P_2P_GS.csv" "Schaef_300_FIX_24P_2P_GS_cens.csv")
counter=0
for proc_file in ${fmri_files[@]}; do
	echo $proc_file
	name_to_parc=${FC_files[$counter]}

	matlab -nodisplay -r "bandpass_filter('$func_loc', '"$proc_file".nii.gz', '$brain_mask', '$TR', '0.008', '0.08', '$repo_directory'); exit;"
	
	fslmaths $func_loc/"$proc_file"_bpass.nii.gz -mul "$gm".nii.gz $func_loc/epi_w
	fslmeants -i $func_loc/epi_w --label=$parc_file -o $func_loc/fsl_parc_ts.txt
	matlab -nodisplay -r "linden_parc('$func_loc', '$name_to_parc'); exit;"
	counter=$((counter+1))
done

rm -r $func_loc/tmp*
#find $func_loc -type f -name "*.nii.gz" ! -name "filtered_func_data.nii.gz" -exec rm {} +
#rm $func_loc/*.nii

