function run_KRR_repeats(FC_name,y_name, covars, outdir,orig_outstem)
%% Run KRR over 20 different seeds
% 
% FC_name = file path to a .mat file where the FC file is saved as roi x
% roi x subjects
%
% y_name = file path to a csv/txt of predictors with headings (should be in the same subject 
% order as the FC)
%
% covars = file path to a csv/txt of covariates you want to include in the
% analysis (no headers, subject order should be the same as FC_name and
% y_name)
% 
% outdir = directory where results will go
%
% orig_outstem= string that results will be saved with in the directory
	
	% Manual Inputs
	repo_directory=; % a string containing the directory the repo is located
    
    data_set='HCP';
    family_csv=('restricted_1110.csv');
    subject_list=('~/kg98_scratch/Kane/KRR/final_KRR_files/HCP/final_sub_list_censrmv.csv');

	%% Start of script
	% Addpaths
    cd repo_directory
    addpath(genpath(sprintf('%s/utils/CBIG-master',repo_directory)))
    
 	data_csv=y_name;
    RSFC_file=(FC_name);
    num_test_folds=20;
    num_inner_folds=20;
    seeds=[5,10,42,33,104,6,27,58,64,82,99,107,113,129,184,172,155,146,111,163];
    orig_outstem=string(orig_outstem);
%% Read predictor variables

        setup_param.y = readtable(data_csv);
        setup_param.y = table2array(setup_param.y(:,2:end));

%% Set up other parameters

        setup_param.covariates = dlmread(covars);
        if data_set=='HCP'
            setup_param.feature_mat = load(RSFC_file); 
            fns=fieldnames(setup_param.feature_mat);
            setup_param.feature_mat = setup_param.feature_mat.(fns{1});
            %FC_idx=dlmread('FC_idx_file); 
            %setup_param.feature_mat=setup_param.feature_mat(:,:,logical(FC_idx));
        elseif data_set=='GSP'
            setup_param.feature_mat = load(RSFC_file); 
            fns=fieldnames(setup_param.feature_mat);
            setup_param.feature_mat = setup_param.feature_mat.(fns{1});
            %FC_idx=dlmread(FC_idx_file); 
            %setup_param.feature_mat=setup_param.feature_mat(:,:,logical(FC_idx));
        end
        setup_param.num_inner_folds = num_inner_folds;
        setup_param.outdir = outdir;
        setup_param.with_bias = 1;
        setup_param.ker_param.type = 'corr';
        setup_param.ker_param.scale = NaN;
        setup_param.lambda_set=load('stable_projects/preprocessing/Li2019_GSR/KernelRidgeRegression/lambda_set.mat');
        setup_param.lambda_set = setup_param.lambda_set.lambda_set;
        setup_param.threshold_set = [];
        setup_param.metric='corr';
        setup_param.cov_X=[];
%% CV split

    for rep=1:length(seeds)
        seed=seeds(rep);
        setup_param.outstem=char(orig_outstem + '_rep' +num2str(rep));
        if data_set=='HCP'
            CBIG_cross_validation_data_split( subject_list, family_csv, 'Subject', ...
                'Family_ID', num_test_folds, seed, fullfile(outdir, ['randseed_' num2str(seed)]), ',' );
        else
            CBIG_cross_validation_data_split( subject_list, 'NONE', 'NONE', ...
            'NONE', num_test_folds, seed, fullfile(outdir, ['randseed_' num2str(seed)]), ',' );
        end
        sub_fold_file = fullfile(outdir, ['randseed_' num2str(seed)], ...
            ['no_relative_' num2str(num_test_folds) '_fold_sub_list.mat']);
        sub_fold = load(sub_fold_file);
        setup_param.sub_fold = sub_fold.sub_fold;


%% Run KRR

        CBIG_KRR_workflow_LITE(setup_param)
    end
end

%% Collate results
% cd ~/kg98_scratch/Kane/KRR/GSP
% procs={'2P', 'AROMA', 'AROMA_2P', 'GSR', 'Dic1', 'Dic2', 'Dic3', 'Dic4', 'Dic5', 'Dic10', 'indiv_dicer'};
% results=struct();
% fnames={'AROMA1','AROMA2','AROMA3','AROMA4','AROMA5','AROMA6','AROMA7','AROMA8','AROMA9','AROMA10'};
% for a=1:10
%     for p=1:length(procs)
%         curr_proc=string(procs{p});
%         if p==1
%             res_name=sprintf('%s/final_result_%s.mat', curr_proc,curr_proc);
%             load(res_name);
%         else
%             res_name=sprintf('AROMA%d/%s/final_result_%s',a,curr_proc,curr_proc);
%             load(res_name);
%         end
%         curr_res(:,:,p)=optimal_acc;
%     end
%     results.(fnames{a})=curr_res;
% end
% save('~/kg98_scratch/Kane/GSP/Results/KRR/all_accuracies.mat','results')


%% For FIX HCP
%cd ~/kg98_scratch/Kane/KRR/HCP/FINAL_FINAL
%procs={"24P", "24P_2P","FIX_24P_2P", "FIX_24P_2P_cens", "FIX_24P_2P_GS", "FIX_24P_2P_GS_cens", "FIX_24P_2P_Dic1", "FIX_24P_2P_Dic1_cens", "FIX_24P_2P_Dic2", "FIX_24P_2P_Dic2_cens", "FIX_24P_2P_Dic3" %"FIX_24P_2P_Dic3_cens", "FIX_24P_2P_Dic4", "FIX_24P_2P_Dic4_cens" "FIX_24P_2P_Dic5", "FIX_24P_2P_Dic5_cens", "FIX_24P_2P_DicOpt", "FIX_24P_2P_DicOpt_cens"};

%for p=1:length(procs)
%    curr_proc=string(procs{p});
%    for rep=1:20
%        res_name=sprintf('%s/final_result_%s_rep%d.mat', curr_proc,curr_proc,rep);
%        load(res_name);
%        curr_rep_res(:,:,rep)=optimal_acc;
%    end

%    results(:,:,:,p)=curr_rep_res;
%end
%save('~/kg98_scratch/Kane/HCP/Results/KRR/all_accuracies.mat','results')

%% For FIX GSP
% cd ~/kg98_scratch/Kane/KRR/GSP
% procs={'2P','2P_GSR', 'FIX', 'FIX_2P', 'GSR', 'Dic1', 'Dic2', 'Dic3', 'Dic4', 'Dic5', 'Dic10', 'indiv_dicer'};
% results=struct();
%     for p=1:length(procs)
%         curr_proc=string(procs{p});
%         if strcmp(curr_proc, '2P') || strcmp(curr_proc, '2P_GSR')
%             res_name=sprintf('%s/final_result_%s.mat', curr_proc,curr_proc);
%             load(res_name);
%         else
%             res_name=sprintf('FIX/%s/final_result_%s.mat', curr_proc,curr_proc);
%             load(res_name);
%         end
%         curr_res(:,:,p)=optimal_acc;
%     end
% results.FIX=curr_res;
% save('~/kg98_scratch/Kane/GSP/FIX/Results/KRR/all_accuracies_w2pgsr.mat','results')
