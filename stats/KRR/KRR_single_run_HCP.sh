#!/bin/bash

#SBATCH --job-name=KRR
#SBATCH --time=0-6:00:00
#SBATCH --ntasks=1
#SBATCH --mem=20000
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mail-user=<knpavlovich@gmail.com>
#SBATCH --mail-type=FAIL
while getopts "s:f:d:" flag; do
	case "${flag}" in
		s) outstem=${OPTARG} ;;
		f) FC_name=${OPTARG} ;;
		d) outdir=${OPTARG} ;;
	esac
done
module purge
module load matlab

cd ~/kg98_scratch/Kane/KRR

y_name="final_KRR_files/HCP/final_58_y_list_censrmv.csv"
covars="final_KRR_files/HCP/final_covars_gend_age_fd_censrmv.csv"

echo $FC_name
matlab -nodisplay -r "run_KRR_repeats('$FC_name','$y_name','$covars','$outdir','$outstem'); exit;"


