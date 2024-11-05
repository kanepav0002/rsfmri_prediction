#!/bin/bash

#SBATCH --job-name=GSP_FIX_dicer
#SBATCH --time=0-8:30:00
#SBATCH --ntasks=1
#SBATCH --mem=65000
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mail-user=<knpavlovich@gmail.com>
#SBATCH --mail-type=FAIL
while getopts "s:r:" flag; do
	case "${flag}" in
		s) sub=${OPTARG} ;; 
		r) reps=${OPTARG} ;;
	esac
done

# Note for this script to work you need to extract the motion parameters from fmriprep and put them in the func folder labelled as prep_motion_params.par
# Also need a list of censored frames that you can get from the script get_censored_frames.m

# Manual inputs
repo_directory='location the repo has been downloaded to'
func_loc="subjects functional directory e.g. /home/CNP/$sub/func"
anat_loc="subjects anatomical directory e.g. /home/CNP/$sub/anat"
TR="TR of the fmri scan in seconds e.g. 3"
TR_ms="TR of the fmri scan in ms e.g. 3000"
num_frames="length of fmri scan i.e. the number of frames in the scan"
# Define your own parcellation and the number of parcels below or use the defaults as listed
parc_file="$repo_directory/dat/Schaefer2018_300Parcels_17Networks_order_FSLMNI152_2mm.nii.gz"
num_parcels=300

# Modules required - make sure these exist
# Replace these with the module names in your cluster
# Note that ICA-FIX uses its own modules further down in the code.
module purge 
module load fsl/6.0.1
module load matlab

#activate python
source $repo_directory/DiCER_env/bin/activate

dicer_folder=$repo_directory/utils/DiCER-master/

cd $repo_directory/pre_processing/utils

# Do baseline processing if AROMA script has already been run for the subject no need to do this again.
sh base_proc.sh -s $sub -f $func_loc -a $anat_loc -p $parc_file -r $repo_directory -n $num_frames

brain_mask=$anat_loc/brain_mask.nii.gz
gm="$anat_loc/c1"$sub"_desc-preproc_T1w_brain"
func_file="filtered_func_data_clean.nii.gz"
base_proc="ifr_func.nii.gz"


##############################################################
# Base Processing Pre-FIX
# Will do 24P regression w/ and w/out 2P regression and save FC
#############################################################

cp $func_loc/ifr_func.nii.gz $func_loc/base_proc.nii.gz 

fslmeants -i $func_loc/base_proc.nii.gz -o $func_loc/WM_pre_fix_ts.txt -m $anat_loc/final_wmmask.nii.gz
fslmeants -i $func_loc/base_proc.nii.gz -o $func_loc/CSF_pre_fix_ts.txt -m $anat_loc/final_csfmask.nii.gz

python get_p24_regs.py -s $sub -f $func_loc -w $func_loc/WM_pre_fix_ts.txt -c $func_loc/CSF_pre_fix_ts.txt -g "not_this_time" -d "not_this_time"

# Create a list with no censored frames
dim_time=$(fslinfo $func_loc/base_proc.nii.gz | awk '/dim4/ {print $2}')
dim_time=${dim_time:0:4}
rm $func_loc/no_censor.txt
for ((i=0; i<dim_time; i++)); do
    echo "1" >> $func_loc/no_censor.txt
done

# Do Regression for 24p 
echo $func_loc/base_proc.nii.gz > $func_loc/tmp_fmri_list.txt
echo $func_loc/24P.nii.gz > $func_loc/tmp_out_list.txt
echo $func_loc/p_24_regs.csv > $func_loc/tmp_regs.txt
echo $func_loc/no_censor.txt > $func_loc/tmp_c_frames_list.txt
matlab -nodisplay -r "CBIG_glm_regress_vol('$func_loc/tmp_fmri_list.txt', '$func_loc/tmp_out_list.txt','$func_loc/tmp_regs.txt', '1', '$func_loc/tmp_c_frames_list.txt', '1', '0', '$repo_directory'); exit;"
# And 24 + 2P
echo $func_loc/base_proc.nii.gz > $func_loc/tmp_fmri_list.txt
echo $func_loc/24_2P.nii.gz > $func_loc/tmp_out_list.txt
echo $func_loc/p_24_2p_regs.csv > $func_loc/tmp_regs.txt
echo $func_loc/no_censor.txt > $func_loc/tmp_c_frames_list.txt
matlab -nodisplay -r "CBIG_glm_regress_vol('$func_loc/tmp_fmri_list.txt', '$func_loc/tmp_out_list.txt','$func_loc/tmp_regs.txt', '1', '$func_loc/tmp_c_frames_list.txt', '1', '0', '$repo_directory'); exit;"

# Bandpass Filter 
matlab -nodisplay -r "bandpass_filter('$func_loc', '24P.nii.gz', '$brain_mask', '$TR', '0.008', '0.08', '$repo_directory'); exit;"
matlab -nodisplay -r "bandpass_filter('$func_loc', '24_2P.nii.gz', '$brain_mask', '$TR', '0.008', '0.08', '$repo_directory'); exit;"

# Smooth before parcellating / FIX
matlab -nodisplay -r "smooth('base_proc.nii.gz', '$func_loc', '$repo_directory'); exit;"
matlab -nodisplay -r "smooth('24P_bpass.nii.gz', '$func_loc', '$repo_directory'); exit;"
matlab -nodisplay -r "smooth('24_2P_bpass.nii.gz', '$func_loc', '$repo_directory'); exit;"

# Make FC from those two files
fslmeants -i "$gm".nii --label=$parc_file -o $func_loc/fsl_gm_prob.txt
fslmaths $func_loc/s24P_bpass.nii.gz -mul "$gm".nii $func_loc/epi_w
fslmeants -i $func_loc/epi_w --label=$parc_file -o $func_loc/fsl_parc_ts.txt
matlab -nodisplay -r "linden_parc('$func_loc', 'Schaef300_24P.csv'); exit;"

fslmaths $func_loc/s24_2P_bpass.nii.gz -mul "$gm".nii $func_loc/epi_w
fslmeants -i $func_loc/epi_w --label=$parc_file -o $func_loc/fsl_parc_ts.txt
matlab -nodisplay -r "linden_parc('$func_loc', 'Schaef300_24P_2P.csv'); exit;"


##############################################################
# Do FIX on 24P data
# do one pipeline with Just FIX
# Another where an additional regression of 24p and 2p.
# Split into scrub and no scrub.
# Another with 24p, 2p, GSR
# Another with 24p 2p, DiCER
#############################################################

# Modules required for FIX - make sure these exist
# Replace these with the module names in your cluster
module purge
export LD_LIBRARY_PATH="/usr/lib64/:$LD_LIBRARY_PATH"
export R_LIBS="$repo_directory/utils/FIX_R_LIB"
module load R/4.4.0-mkl
module load matlab/r2017a
module load fsl/5.0.9
# FIX the data
cd $func_loc
cd ../../ # Fix requires you to be in directory above subjects directory

sh $repo_directory/pre_processing/utils/prep_FIX_files.sh -s $sub -f $func_loc -r $repo_directory
sh $repo_directory/utils/fix/fix $sub $repo_directory/dat/GSP_Train_Weights.RData 40

cd $repo_directory/pre_processing/utils

cp $func_loc/../filtered_func_data_clean.nii.gz $func_loc/filtered_func_data_clean.nii.gz
# Modules required for the rest
module purge 
module load fsl/6.0.1
source /home/kanep/kg98_scratch/Kane/DiCERxModularity/DiCER_env/bin/activate
module load matlab

# Get new regs
fslmeants -i $func_loc/filtered_func_data_clean.nii.gz -o $func_loc/WM_post_fix_ts.txt -m $anat_loc/final_wmmask.nii.gz
fslmeants -i $func_loc/filtered_func_data_clean.nii.gz -o $func_loc/CSF_post_fix_ts.txt -m $anat_loc/final_csfmask.nii.gz
fslmeants -i $func_loc/filtered_func_data_clean.nii.gz -o $func_loc/GS_post_fix_ts.txt -m $anat_loc/brain_mask.nii.gz

python get_p24_regs.py -s $sub -f $func_loc -w $func_loc/WM_post_fix_ts.txt -c $func_loc/CSF_post_fix_ts.txt -g $func_loc/GS_post_fix_ts.txt -d "not_this_time"

# Do Regressions with and without ignoring scrubbed frames
echo $func_loc/filtered_func_data_clean.nii.gz > $func_loc/tmp_fmri_list.txt
echo $func_loc/FIX_24P_2P.nii.gz > $func_loc/tmp_out_list.txt
echo $func_loc/p_24_2p_regs.csv > $func_loc/tmp_regs.txt
echo $func_loc/no_censor.txt > $func_loc/tmp_c_frames_list.txt
matlab -nodisplay -r "CBIG_glm_regress_vol('$func_loc/tmp_fmri_list.txt', '$func_loc/tmp_out_list.txt','$func_loc/tmp_regs.txt', '1', '$func_loc/tmp_c_frames_list.txt', '1', '0', '$repo_directory'); exit;"
echo $func_loc/FIX_24P_2P_c.nii.gz > $func_loc/tmp_out_list.txt
echo $func_loc/censored_frames.txt > $func_loc/tmp_c_frames_list.txt
matlab -nodisplay -r "CBIG_glm_regress_vol('$func_loc/tmp_fmri_list.txt', '$func_loc/tmp_out_list.txt','$func_loc/tmp_regs.txt', '1', '$func_loc/tmp_c_frames_list.txt', '1', '0', '$repo_directory'); exit;"

echo $func_loc/filtered_func_data_clean.nii.gz > $func_loc/tmp_fmri_list.txt
echo $func_loc/FIX_24P_2P_GS.nii.gz > $func_loc/tmp_out_list.txt
echo $func_loc/p_24_2p_GS_regs.csv > $func_loc/tmp_regs.txt
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
sh DiCER_very_lightweight.sh -i FIX_24P_2P.nii.gz  -t ../../$sub/anat/"$sub"_dtissue_func.nii.gz -w $func_loc -s $sub -r $reps -e 0.8 

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
	
	fslmaths $func_loc/FIX_24P_2P_Dic"$dic_num"_bpass.nii.gz -mul "$gm".nii $func_loc/epi_w
	fslmeants -i $func_loc/epi_w --label=$parc_file -o $func_loc/fsl_parc_ts.txt
	matlab -nodisplay -r "linden_parc('$func_loc', 'Schaef_300_FIX_24P_2P_Dic"$dic_num".csv'); exit;"

	fslmaths $func_loc/FIX_24P_2P_Dic"$dic_num"_cens_bpass.nii.gz -mul "$gm".nii $func_loc/epi_w
	fslmeants -i $func_loc/epi_w --label=$parc_file -o $func_loc/fsl_parc_ts.txt
	matlab -nodisplay -r "linden_parc('$func_loc', 'Schaef_300_FIX_24P_2P_Dic"$dic_num"_cens.csv'); exit;"
done

# Get modularity Maximised FC
python mod_max_FC.py -f $func_loc -r $reps -b "FC_Schaef_300_FIX_24P_2P_Dic" -p $num_parcels

# bandpass filter and parcellate all FIX files
fmri_files=("filtered_func_data" "FIX_24P_2P" "FIX_24P_2P_cens" "FIX_24P_2P_GS" "FIX_24P_2p_GS_cens")
FC_files=("Schaef300_FIX.csv" "Schaef_300_FIX_24P_2P.csv" "Schaef_300_FIX_24P_2P_cens.csv" "Schaef_300_FIX_24P_2P_GS.csv" "Schaef_300_FIX_24P_2p_GS_cens.csv")
counter=0
for proc_file in ${fmri_files[@]}; do
	echo $proc_file
	name_to_parc=${FC_files[$counter]}

	matlab -nodisplay -r "bandpass_filter('$func_loc', '"$proc_file".nii.gz', '$brain_mask', '$TR', '0.008', '0.08', '$repo_directory'); exit;"
	
	fslmaths $func_loc/"$proc_file"_bpass.nii.gz -mul "$gm".nii $func_loc/epi_w
	fslmeants -i $func_loc/epi_w --label=$parc_file -o $func_loc/fsl_parc_ts.txt
	matlab -nodisplay -r "linden_parc('$func_loc', '$name_to_parc'); exit;"
	counter=$((counter+1))
done


rm -r $func_loc/tmp*
#rm $func_loc/*.nii
#find $func_loc -type f -name "*.nii.gz" ! -name "filtered_func_data_clean.nii.gz" -exec rm {} +

