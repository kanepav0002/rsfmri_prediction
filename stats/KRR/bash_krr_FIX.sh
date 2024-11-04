#!/bin/bash

#SBATCH --job-name=GSP_KRR
#SBATCH --time=0-12:30:00
#SBATCH --ntasks=1
#SBATCH --mem=45000
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mail-user=<knpavlovich@gmail.com>
#SBATCH --mail-type=FAIL

module purge
module load matlab

cd ~/kg98_scratch/Kane/KRR

FC_files=("../HCP/Results/final_final_final_FC/24P" "../HCP/Results/final_final_final_FC/24P_2P" "../HCP/Results/final_final_final_FC/FIX" "../HCP/Results/final_final_final_FC/FIX_24P_2P" "../HCP/Results/final_final_final_FC/FIX_24P_2P_cens" "../HCP/Results/final_final_final_FC/FIX_24P_2P_GS" "../HCP/Results/final_final_final_FC/FIX_24P_2P_GS_cens" "../HCP/Results/final_final_final_FC/FIX_24P_2P_Dic1" "../HCP/Results/final_final_final_FC/FIX_24P_2P_Dic1_cens" "../HCP/Results/final_final_final_FC/FIX_24P_2P_Dic2" "../HCP/Results/final_final_final_FC/FIX_24P_2P_Dic2_cens" "../HCP/Results/final_final_final_FC/FIX_24P_2P_Dic3.mat" "../HCP/Results/final_final_final_FC/FIX_24P_2P_Dic3_cens" "../HCP/Results/final_final_final_FC/FIX_24P_2P_Dic4" "../HCP/Results/final_final_final_FC/FIX_24P_2P_Dic4_cens" "../HCP/Results/final_final_final_FC/FIX_24P_2P_Dic5" "../HCP/Results/final_final_final_FC/FIX_24P_2P_Dic5_cens" "../HCP/Results/final_final_final_FC/FIX_24P_2P_DicOpt" "../HCP/Results/final_final_final_FC/FIX_24P_2P_DicOpt_cens")

stem_names=("24P" "24P_2P" "FIX" "FIX_24P_2P" "FIX_24P_2P_cens" "FIX_24P_2P_GS" "FIX_24P_2P_GS_cens" "FIX_24P_2P_Dic1" "FIX_24P_2P_Dic1_cens" "FIX_24P_2P_Dic2" "FIX_24P_2P_Dic2_cens" "FIX_24P_2P_Dic3" "FIX_24P_2P_Dic3_cens" "FIX_24P_2P_Dic4" "FIX_24P_2P_Dic4_cens" "FIX_24P_2P_Dic5" "FIX_24P_2P_Dic5_cens" "FIX_24P_2P_DicOpt" "FIX_24P_2P_DicOpt_cens")

counter=0
for FC in ${FC_files[@]}; do
	
	FC_name=${FC_files[$counter]}
	outdir="HCP/FINAL_FINAL/"${stem_names[$counter]}
	outstem=${stem_names[$counter]}
	counter=$((counter+1))

	sbatch KRR_single_run_HCP.sh -s $outstem -f $FC_name -d $outdir
done

