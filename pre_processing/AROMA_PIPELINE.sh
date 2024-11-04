#!/bin/bash

#SBATCH --job-name=AROMA Pipeline
#SBATCH --time=0-24:30:00
#SBATCH --ntasks=1
#SBATCH --mem=65000
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

while getopts "s:r:" flag; do
	case "${flag}" in
		s) sub=${OPTARG} ;; 
		r) reps=${OPTARG} ;;
	esac
done

# Note for this script to work you need to extract the motion parameters from fmriprep and put them in the func folder labelled as prep_motion_params.par
# Also need a list of censored frames
# This can be done with the scripts in the folder pre-pre-processing

# Manual inputs
repo_directory='location the repo has been downloaded to'
func_loc="subjects functional directory e.g. /home/CNP/$sub/func"
anat_loc="subjects anatomical directory e.g. /home/CNP/$sub/anat"
TR="TR of the fmri scan in seconds e.g. 3"
TR_ms="TR of the fmri scan in ms e.g. 3000"
num_frames="number of frames in the fmri file e.g. 152"
# Below can be changed if you want to use your own parcellation, otherwise leave defaults below
parc_file="$repo_directory/dat/Schaefer2018_300Parcels_17Networks_order_FSLMNI152_2mm.nii.gz"
num_parcels=300

# Modules required - make sure these exist
# Replace these with the module names in your cluster
module purge 
module load fsl/6.0.1
module load matlab

#activate python
source $repo_directory/DiCER_env/bin/activate

dicer_folder=$repo_directory/utils/DiCER-master/

cd $repo_directory/pre_processing/utils

# Do baseline processing if FIX script has already been run for GSP no need to do this again.
sh base_proc.sh -s $sub -f $func_loc -a $anat_loc -p $parc_file -r $repo_directory -n $num_frames

brain_mask=$anat_loc/brain_mask.nii.gz
gm="$anat_loc/c1"$sub"_desc-preproc_T1w_brain"
base_proc="ifr_func.nii.gz"

##############################################################
# Base Processing Pre-AROMA
# Will do 24P regression w/ and w/out 2P regression and save FC
#############################################################

fslmeants -i $func_loc/$base_proc -o $func_loc/WM_pre_fix_ts.txt -m $anat_loc/final_wmmask.nii.gz
fslmeants -i $func_loc/$base_proc -o $func_loc/CSF_pre_fix_ts.txt -m $anat_loc/final_csfmask.nii.gz

python get_p24_regs.py -s $sub -f $func_loc -w $func_loc/WM_pre_fix_ts.txt -c $func_loc/CSF_pre_fix_ts.txt -g "not_this_time" -d "not_this_time"

# Create a list with no censored frames
dim_time=$(fslinfo $func_loc/$base_proc | awk '/dim4/ {print $2}')
dim_time=${dim_time:0:4}
rm $func_loc/no_censor.txt
for ((i=0; i<dim_time; i++)); do
    echo "1" >> $func_loc/no_censor.txt
done

# Do Regression for 24p 
echo $func_loc/$base_proc > $func_loc/tmp_fmri_list.txt
echo $func_loc/24P.nii.gz > $func_loc/tmp_out_list.txt
echo $func_loc/p_24_regs.csv > $func_loc/tmp_regs.txt
echo $func_loc/no_censor.txt > $func_loc/tmp_c_frames_list.txt
matlab -nodisplay -r "CBIG_glm_regress_vol('$func_loc/tmp_fmri_list.txt', '$func_loc/tmp_out_list.txt','$func_loc/tmp_regs.txt', '1', '$func_loc/tmp_c_frames_list.txt', '1', '0', '$repo_directory'); exit;"
# And 24 + 2P
echo $func_loc/$base_proc > $func_loc/tmp_fmri_list.txt
echo $func_loc/24P_2P.nii.gz > $func_loc/tmp_out_list.txt
echo $func_loc/p_24_2p_regs.csv > $func_loc/tmp_regs.txt
echo $func_loc/no_censor.txt > $func_loc/tmp_c_frames_list.txt
matlab -nodisplay -r "CBIG_glm_regress_vol('$func_loc/tmp_fmri_list.txt', '$func_loc/tmp_out_list.txt','$func_loc/tmp_regs.txt', '1', '$func_loc/tmp_c_frames_list.txt', '1', '0', '$repo_directory'); exit;"

# Bandpass Filter 
matlab -nodisplay -r "bandpass_filter('$func_loc', '24P.nii.gz', '$brain_mask', '$TR', '0.008', '0.08', '$repo_directory'); exit;"
matlab -nodisplay -r "bandpass_filter('$func_loc', '24P_2P.nii.gz', '$brain_mask', '$TR', '0.008', '0.08', '$repo_directory'); exit;"

# Smooth before parcellating / AROMA
matlab -nodisplay -r "smooth('$base_proc', '$func_loc', '$repo_directory'); exit;"
matlab -nodisplay -r "smooth('24P_bpass.nii.gz', '$func_loc',  '$repo_directory'); exit;"
matlab -nodisplay -r "smooth('24P_2P_bpass.nii.gz', '$func_loc',  '$repo_directory'); exit;"

# Make FC from those two files
fslmeants -i "$gm".nii --label=$parc_file -o $func_loc/fsl_gm_prob.txt
fslmaths $func_loc/s24P_bpass.nii.gz -mul "$gm".nii $func_loc/epi_w
fslmeants -i $func_loc/epi_w --label=$parc_file -o $func_loc/fsl_parc_ts.txt
matlab -nodisplay -r "linden_parc('$func_loc', 'Schaef300_24P.csv'); exit;"

fslmaths $func_loc/s24P_2P_bpass.nii.gz -mul "$gm".nii $func_loc/epi_w
fslmeants -i $func_loc/epi_w --label=$parc_file -o $func_loc/fsl_parc_ts.txt
matlab -nodisplay -r "linden_parc('$func_loc', 'Schaef300_24P_2P.csv'); exit;"

# Remove repeat files
rm $func_loc/*.nii
rm $func_loc/*.mat
##############################################################
# Do AROMA with and without 24P 
# Another where an additional regression of 2p.
# Split into scrub and no scrub.
# Another with AROMA, 24p, 2p, GSR
# Another with AROMA, 24p 2p, 1DiCER
# Repeat all AROMA pipelines 10 times to assess variation across runs
#############################################################

for a in $(seq 1 2)
do
	AROMA_dir="AROMA"$a""
	python $repo_directory/utils/ICA-AROMA-master/ICA_AROMA.py -in $func_loc/"s"$base_proc -out $func_loc/$AROMA_dir/ -mc $func_loc/prep_motion_params.par -m $brain_mask -tr $TR -den nonaggr -overwrite
	
	# Get new regs
	fslmeants -i $func_loc/$AROMA_dir/denoised_func_data_nonaggr.nii.gz -o $func_loc/$AROMA_dir/WM_post_AROMA_ts.txt -m $anat_loc/final_wmmask.nii.gz
	fslmeants -i $func_loc/$AROMA_dir/denoised_func_data_nonaggr.nii.gz -o $func_loc/$AROMA_dir/CSF_post_AROMA_ts.txt -m $anat_loc/final_csfmask.nii.gz
	fslmeants -i $func_loc/$AROMA_dir/denoised_func_data_nonaggr.nii.gz -o $func_loc/$AROMA_dir/GS_post_AROMA_ts.txt -m $anat_loc/brain_mask.nii.gz

	# Collate Regressors and get first derivatives.
	cp $func_loc/prep_motion_params.par $func_loc/$AROMA_dir/
	cp $func_loc/fsl_gm_prob.txt $func_loc/$AROMA_dir/fsl_gm_prob.txt
	python get_p24_regs.py -s $sub -f $func_loc/$AROMA_dir -w $func_loc/$AROMA_dir/WM_post_AROMA_ts.txt -c $func_loc/$AROMA_dir/CSF_post_AROMA_ts.txt -g $func_loc/$AROMA_dir/GS_post_AROMA_ts.txt -d $func_loc/$AROMA_dir/"$sub"_dbscan_1regressors.tsv

	# Do Regressions with and without ignoring scrubbed frames
	echo $func_loc/$AROMA_dir/denoised_func_data_nonaggr.nii.gz > $func_loc/tmp_fmri_list.txt
	echo $func_loc/$AROMA_dir/AROMA_24P_2P.nii.gz > $func_loc/tmp_out_list.txt
	echo $func_loc/$AROMA_dir/p_24_2p_regs.csv > $func_loc/tmp_regs.txt
	echo $func_loc/no_censor.txt > $func_loc/tmp_c_frames_list.txt
	matlab -nodisplay -r "CBIG_glm_regress_vol('$func_loc/tmp_fmri_list.txt', '$func_loc/tmp_out_list.txt','$func_loc/tmp_regs.txt', '1', '$func_loc/tmp_c_frames_list.txt','1', '0', '$repo_directory'); exit;"
	echo $func_loc/$AROMA_dir/AROMA_24P_2P_c.nii.gz > $func_loc/tmp_out_list.txt
	echo $func_loc/censored_frames.txt > $func_loc/tmp_c_frames_list.txt
	matlab -nodisplay -r "CBIG_glm_regress_vol('$func_loc/tmp_fmri_list.txt', '$func_loc/tmp_out_list.txt','$func_loc/tmp_regs.txt', '1', '$func_loc/tmp_c_frames_list.txt', '1', '0', '$repo_directory'); exit;"

	echo $func_loc/$AROMA_dir/denoised_func_data_nonaggr.nii.gz > $func_loc/tmp_fmri_list.txt
	echo $func_loc/$AROMA_dir/AROMA_24P_2P_GS.nii.gz > $func_loc/tmp_out_list.txt
	echo $func_loc/$AROMA_dir/p_24_2p_GS_regs.csv > $func_loc/tmp_regs.txt
	echo $func_loc/no_censor.txt > $func_loc/tmp_c_frames_list.txt
	matlab -nodisplay -r "CBIG_glm_regress_vol('$func_loc/tmp_fmri_list.txt', '$func_loc/tmp_out_list.txt','$func_loc/tmp_regs.txt', '1', '$func_loc/tmp_c_frames_list.txt','1', '0', '$repo_directory'); exit;"
	echo $func_loc/$AROMA_dir/AROMA_24P_2P_GS_c.nii.gz > $func_loc/tmp_out_list.txt
	echo $func_loc/censored_frames.txt > $func_loc/tmp_c_frames_list.txt
	matlab -nodisplay -r "CBIG_glm_regress_vol('$func_loc/tmp_fmri_list.txt', '$func_loc/tmp_out_list.txt','$func_loc/tmp_regs.txt', '1', '$func_loc/tmp_c_frames_list.txt','1', '0', '$repo_directory'); exit;"

	# Scrub the data
	matlab -nodisplay -r "CBIG_preproc_censor_wrapper('$func_loc/$AROMA_dir/AROMA_24P_2P_GS_c.nii.gz', '$func_loc/censored_frames.txt', '$TR_ms', '$func_loc/tmp_interp_bold.nii.gz', '$func_loc/$AROMA_dir/AROMA_24P_2p_GS_cens.nii.gz', '$brain_mask', '10', '$repo_directory'); exit;"
	matlab -nodisplay -r "CBIG_preproc_censor_wrapper('$func_loc/$AROMA_dir/AROMA_24P_2P_c.nii.gz', '$func_loc/censored_frames.txt', '$TR_ms', '$func_loc/tmp_interp_bold.nii.gz', '$func_loc/$AROMA_dir/AROMA_24P_2P_cens.nii.gz', '$brain_mask', '10', '$repo_directory'); exit;"
	
# Do DiCER
	cd $dicer_folder
	echo "Starting DiCER"
	sh DiCER_very_lightweight.sh -i AROMA_24P_2P.nii.gz  -t ../../anat/"$sub"_dtissue_func.nii.gz -w $func_loc/$AROMA_dir -s $sub -r $reps -e 0.8 -d
	cd $repo_directory/pre_processing/utils
	for dic_num in $(seq 1 $reps); do

		python get_dic_regs.py -s $sub -f $func_loc/$AROMA_dir -d $func_loc/$AROMA_dir/"$sub"_dbscan_"$dic_num"regressors.tsv

		echo $func_loc/$AROMA_dir/AROMA_24P_2P.nii.gz > $func_loc/tmp_fmri_list.txt
		echo $func_loc/$AROMA_dir/AROMA_24P_2P_Dic"$dic_num".nii.gz > $func_loc/tmp_out_list.txt
		echo $func_loc/$AROMA_dir/Dic_regs.csv > $func_loc/tmp_regs.txt
		echo $func_loc/no_censor.txt > $func_loc/tmp_c_frames_list.txt
		matlab -nodisplay -r "CBIG_glm_regress_vol('$func_loc/tmp_fmri_list.txt', '$func_loc/tmp_out_list.txt','$func_loc/tmp_regs.txt', '1', '$func_loc/tmp_c_frames_list.txt','1', '0', '$repo_directory'); exit;"
		echo $func_loc/$AROMA_dir/AROMA_24P_2P_Dic"$dic_num"_c.nii.gz > $func_loc/tmp_out_list.txt
		echo $func_loc/censored_frames.txt > $func_loc/tmp_c_frames_list.txt
		matlab -nodisplay -r "CBIG_glm_regress_vol('$func_loc/tmp_fmri_list.txt', '$func_loc/tmp_out_list.txt','$func_loc/tmp_regs.txt', '1', '$func_loc/tmp_c_frames_list.txt','1', '0', '$repo_directory'); exit;"

		matlab -nodisplay -r "CBIG_preproc_censor_wrapper('$func_loc/$AROMA_dir/AROMA_24P_2P_Dic"$dic_num"_c.nii.gz', '$func_loc/censored_frames.txt', '$TR_ms', '$func_loc/tmp_interp_bold.nii.gz', '$func_loc/$AROMA_dir/AROMA_24P_2P_Dic"$dic_num"_cens.nii.gz', '$brain_mask', '10', '$repo_directory'); exit;"

		matlab -nodisplay -r "bandpass_filter('$func_loc/$AROMA_dir', 'AROMA_24P_2P_Dic"$dic_num".nii.gz', '$brain_mask', '$TR', '0.008', '0.08', '$repo_directory'); exit;"
		matlab -nodisplay -r "bandpass_filter('$func_loc/$AROMA_dir', 'AROMA_24P_2P_Dic"$dic_num"_cens.nii.gz', '$brain_mask', '$TR', '0.008', '0.08', '$repo_directory'); exit;"
	
		fslmaths $func_loc/$AROMA_dir/AROMA_24P_2P_Dic"$dic_num"_bpass.nii.gz -mul "$gm".nii $func_loc/$AROMA_dir/epi_w
		fslmeants -i $func_loc/$AROMA_dir/epi_w --label=$parc_file -o $func_loc/$AROMA_dir/fsl_parc_ts.txt
		matlab -nodisplay -r "linden_parc('$func_loc/$AROMA_dir', 'Schaef_300_AROMA_24P_2P_Dic"$dic_num".csv'); exit;"

		fslmaths $func_loc/$AROMA_dir/AROMA_24P_2P_Dic"$dic_num"_cens_bpass.nii.gz -mul "$gm".nii $func_loc/$AROMA_dir/epi_w
		fslmeants -i $func_loc/$AROMA_dir/epi_w --label=$parc_file -o $func_loc/$AROMA_dir/fsl_parc_ts.txt
		matlab -nodisplay -r "linden_parc('$func_loc/$AROMA_dir', 'Schaef_300_AROMA_24P_2P_Dic"$dic_num"_cens.csv'); exit;"
	done
	
	# Get modularity Maximised FC
	python mod_max_FC.py -f $func_loc/$AROMA_dir -r $reps -b "FC_Schaef_300_AROMA_24P_2P_Dic" -p $num_parcels
	# bandpass filter and parcellate all FIX files
	fmri_files=("denoised_func_data_nonaggr" "AROMA_24P_2P" "AROMA_24P_2P_cens" "AROMA_24P_2P_GS" "AROMA_24P_2p_GS_cens")
	FC_files=("Schaef300_AROMA.csv" "Schaef300_AROMA_24P_2P.csv" "Schaef300_AROMA_24P_2P_cens.csv" "Schaef300_AROMA_24P_2P_GS.csv" "Schaef300_AROMA_24P_2p_GS_cens.csv")
	counter=0
	for proc_file in ${fmri_files[@]}; do
		echo $proc_file
		name_to_parc=${FC_files[$counter]}

		matlab -nodisplay -r "bandpass_filter('$func_loc/$AROMA_dir', '"$proc_file".nii.gz', '$brain_mask', '$TR', '0.008', '0.08', '$repo_directory'); exit;"
	
		fslmeants -i "$gm".nii --label=$parc_file -o $func_loc/$AROMA_dir/fsl_gm_prob.txt
		fslmaths $func_loc/$AROMA_dir/"$proc_file"_bpass.nii.gz -mul "$gm".nii $func_loc/$AROMA_dir/epi_w
		fslmeants -i $func_loc/$AROMA_dir/epi_w --label=$parc_file -o $func_loc/$AROMA_dir/fsl_parc_ts.txt
		matlab -nodisplay -r "linden_parc('$func_loc/$AROMA_dir', '$name_to_parc'); exit;"
		counter=$((counter+1))
	done

	rm -r $func_loc/$AROMA_dir/tmp*
	#rm $func_loc/$AROMA_dir/*.nii*
done
rm -r $func_loc/tmp*
#find $func_loc -type f -name "*.nii.gz" ! -name "$base_proc" -exec rm {} +

