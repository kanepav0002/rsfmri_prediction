This repository replicates the pre-processing and analysis that is found in: 
Pavlovich, K., Pang, J., Holmes, A., Constable, T., & Fornito, A. (2024) The efficacy of resting-state fMRI denoising pipelines for motion correction and behavioural prediction. 

## Requirements
All analysis has been run using versions of the following software
FSL (6.0.1)
MATLAB (r2023b)
R (4.4.0)

## Pre-pre-processing
For CNP and GSP data fmriprep needs to be run on the data before applying the pipeline scripts. 
The scripts will work with the outputs produced by fmriprep, the files needed are:
$sub/func/""$sub"_task-rest_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz"
$sub/anat/"$sub"_desc-preproc_T1w.nii.gz"
$sub/func/""$sub"_task-rest_desc-confounds_timeseries.tsv"

For HCP the data needs to be organised into a BIDS directory with the nifti rest sessions and T1w image from the MNINonLinear HCP folder.
with the following labels

$sub/func/rfMRI_REST1_LR.nii.gz
$sub/func/rfMRI_REST1_LR_hp2000_clean.nii.gz
$sub/func/Movement_Regressors.txt
$sub/anat/T1w_brain.nii.gz - the anatomical image should be in the same MNI space as the functional images

After running fmriprep, and before the pipelines are run, the motion parameters and framewise displacements need to be extracted for each participant fmripreps list of confounds 
put into each subjects functional directory and labelled as  prep_motion_params.par (for GSP and CNP) an example python script for this is in the utils folder
and is labelled get_motion_from_fprep.py

For HCP the motion parameters are labelled as Movement_Regressors.txt and can just be copied into each subjects func folder.
Framewise displacement adjusted for multiband data can then be calculated with the function mband_powerFD.m

The script get_list_of_censored_frames.m can then be run to return a vector of frames that need to be censored for each subject based on FD thresholds


## Pre-processing
For the ICA-AROMA pipeline there is a standard script called AROMA_PIPELINE.sh that will run the AROMA pipeline on a single subject
For the ICA-FIX Pipelines the pipelines are seperated for the GSP and HCP datasets (GSP_FIX_PIPELINE.sh and HCP_PIPELINE.sh)
as the HCP data has undergone FIX seperately, and the GSP requires training based on weights in the /dat folder.

Each script requires manual input of the dataset location and functional parameters, these can be found at the top of each script.
These scripts can be executed straight from the terminal with the inputs -s (for subject ID) 
and -r (for how many DiCER repitions you wish to do) an example of this is:<br>
`./AROMA_PIPELINE.sh -s "sub-1021" -r 10`

Each subject requires a long time to run, especially for pipelines using ICA-AROMA. If you have access to a cluster, example slurm parameters can be found at the top of each script.

# Stats
## QC_FC and VE1
Once each subject has finished running you can collate the FC from pipelines into a .mat file with the following fields: ROIxROIxsubject
Then you can run the QC_FC_VE1.m script to get all the stats, manual fields are required to be inputted at the top of the script.
you need to collate the iterations required to achieve modularity maximisation into one file. In each subjects directory this is saved in a file $sub/func/opt_iters.txt

## Kernel Ridge Regression
The main script for this can be found in stats/KRR/run_KRR_repeats.m
This will run a kernel ridge regression 20 seperate times based on the CBIG script CBIG_KRR_workflow_LITE.m a full breakdown of what files it expects and the parameters can be found in this file. 
This kernel ridge regression is courtesy of the CBIG lab, and full credit goes to them! 

I have included an example of how to send this function to a cluster for different processing pipelines in the script bash_KRR_FIX.sh
To collate the accuracies across repitions after this has been run an example can be found at the bottom of the script.


## Tools required
This work is built off the back bone of many scientists before, and not possible without their contributions. For these scripts to work you need to download these repositories and label them
as below into the utils folder of this repository.

These scripts use a lightweight version of DiCER that doesn't reproduce the plots found in the original paper
If you are using DiCER in your own work please cite the paper from the original repository located here: https://github.com/NSBLab/DiCER

I have included a schaefer parcellation for ease of reproducability, but if you are using for your own work please download 
the original parcellation files here: https://github.com/ThomasYeoLab/CBIG/tree/master/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal

CBIG-master: https://github.com/ThomasYeoLab/CBIG<br>
fieldtrip-master: https://github.com/fieldtrip/fieldtrip<br>
fix: https://fsl.fmrib.ox.ac.uk/fsl/docs/#/resting_state/fix_matlab<br>
ICA-AROMA-master: https://github.com/maartenmennes/ICA-AROMA<br>
linden_rsfmri: https://github.com/lindenmp/rs-fMRI<br>
REST-master: https://github.com/Chaogan-Yan/REST<br>
spm8: https://github.com/spm/spm8<br>






