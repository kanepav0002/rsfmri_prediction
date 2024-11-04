cd ~/kg98_scratch/Kane/GSP/FIX
addpath(genpath('~/kg98_scratch/Kane/DiCERxModularity/QC_FC'));
addpath(genpath('~/kg98_scratch/Kane/utils'));

%% Manual Inputs
repo_directory('string of the directory in which the github repository is located');
dataset_directory('string of where the dataset is stored, each subject should have a folder in this directory')
indiv_optimal_k = ; % This is a txt file that has the optimal iterations identified by mod_max.py, will be saved as opt_iters.txt in each subjects file
sub_list = readcell('txt file containing list of subject');
fd_file_name= ; % Name that the Framewise Displacement file is saved as as a string
num_procs=; % Put in the number of pre processing pipelines you are loading
preprocs = ; % Pre processing pipeline names e.g. 
% ["24P", "24P_2P", "FIX", "FIX_24P_2P", "FIX_24P_2P_cens", "FIX_24P_2P_GS", "FIX_24P_2P_GS_cens"... 
%    ,"FIX_24P_2P_Dic1", "FIX_24P_2P_Dic1_cens"...
%    ,"FIX_24P_2P_Dic2", "FIX_24P_2P_Dic2_cens"...
%    ,"FIX_24P_2P_Dic3", "FIX_24P_2P_Dic3_cens"...
%    ,"FIX_24P_2P_Dic4", "FIX_24P_2P_Dic4_cens"...
%    ,"FIX_24P_2P_Dic5", "FIX_24P_2P_Dic5_cens"...
%    ,"FIX_24P_2P_DicOpt", "FIX_24P_2P_DicOpt_cens"]


% Load FC files Example below
%for p=1:length(preprocs)
%    curr_FC=load(sprintf('Results/final_final_final_FC/%s.mat', preprocs(p)));
%    fns=fieldnames(curr_FC);
%    All_FC(:,:,:,p)=curr_FC.(fns{1});
%end


%% Start of script
% add necessary scripts to path
cd dataset_directory
stat_dir=sprintf('%s/stats',repo_directory);
addpath(genpath(stat_dir))

% Calculate FD exclusions
sub_logical_idx=true(size(sub_list));
preprocs = preprocs';
suprathreshold_exclude = false(1,length(sub_list));
over_5mm = false(1,length(sub_list));
over_cens= false(1,length(sub_list));

for i=1:length(sub_list)
    sub=string(sub_list(i));
    
    sub_orig_FD=dlmread(sprintf('../%s/%s',sub,fd_file_name));
    sub_orig_FD(1)=0;
    orig_FDm(i) = (mean(sub_orig_FD));
    
    high_frames=sub_orig_FD>0.3;
    high_frames_percent = sum(high_frames==1) / numel(high_frames);
    if high_frames_percent>0.2
        suprathreshold_exclude(:,i)=true;
    end
    clear high_frames
    if any(sub_orig_FD > 5)
        over_5mm(:,i)=true;
    end
    
    cens_frames=dlmread(sprintf('%s/func/censored_frames.txt', sub));
    n_frames=size(cens_frames,1)/2;
    n_cens=sum(cens_frames==0);
    if n_cens>n_frames
        over_cens(i)=true;
    end
        
end
 

% Exclude subjects with FD > 0.3, any FD > 5, more than 20% > 0.3
over_3mm_thresh = orig_FDm>0.3;
over_FD_thresh = or(over_3mm_thresh, over_5mm);
over_FD_thresh = or(over_FD_thresh, suprathreshold_exclude);
over_FD_thresh = or(over_FD_thresh, over_cens);
under_FD_thresh = ~over_FD_thresh;
orig_FDm = orig_FDm(under_FD_thresh);

All_FC=All_FC(:,:,under_FD_thresh,:);

% Remove any subjects with NaNs
rmv=squeeze(any(any(any(isnan(All_FC),1),2),4));
orig_FDm(rmv)=[];
All_FC(:,:,rmv,:)=[];

All_QC = struct('preprocs',preprocs,...
    'NaNFilter',[],...
    'QCFC',[],...
    'QCFC_corrected_P',[],...
    'prop_sig_corr',[],...
    'prop_sig_corr_corrected',[],...
    'DistDep',[],...
    'DistDep_P',[]);
    
%% Perform QC_FC

for i=1:length(preprocs)
    [QCFC,P] = GetDistCorr(orig_FDm',All_FC(:,:,:,i));
    P = LP_FlatMat(P);
    
    P_corrected = mafdr(P, 'BHFDR','true');
    All_QC(i).QCFC_corrected_P = P_corrected;
    All_QC(i).QCFC = LP_FlatMat(QCFC);
    
    
    All_QC(i).NaNFilter = ~isnan(All_QC(i).QCFC);
    if ~any(All_QC(i).NaNFilter)
        error('FATAL: No data left after filtering NaNs!');
    elseif any(All_QC(i).NaNFilter)
        fprintf(1, '\tRemoved %u NaN samples from data \n',sum(~All_QC(i).NaNFilter));
        All_QC(i).QCFC = All_QC(i).QCFC(All_QC(i).NaNFilter);
        P = P(All_QC(i).NaNFilter);
    end
    
    All_QC(i).prop_sig_corr=round(sum(P<0.05) / numel(P) * 100,2);
    All_QC(i).prop_sig_corr_corrected = round(sum(P_corrected<0.05) / numel(P_corrected) * 100,2);
end

%% Perform Distance Dependant Correlations

ROI_Coords = dlmread(sprintf('%s/stats/roi_Schaef300_MMP.txt',repo_directory));
% the above ROI coords are for Schaefer 300 parcellation, if you are using
% another parcellation you need to generate your own.
ROIDist = pdist2(ROI_Coords,ROI_Coords,'euclidean');
ROIDistVec = LP_FlatMat(ROIDist);

for i=1:length(preprocs)
    [All_QC(i).DistDep, All_QC(i).DistDep_P] = corr(ROIDistVec(All_QC(i).NaNFilter),All_QC(i).QCFC,'type','Spearman');
    
end


  %% Calculate Variance explained by First PC

process_names={'2P', 'FIX', 'FIX_2P', 'FIX_GSR', 'FIX_2P_1_dicer', 'FIX_2P_2_dicer', 'FIX_2P_3_dicer', ... 
    'FIX_2P_4_dicer', 'FIX_2P_5_dicer', 'FIX_2P_10_dicer'};

sub_list=sub_list(under_FD_thresh);
sub_list(rmv)=[];
indiv_optimal_k=indiv_optimal_k(under_FD_thresh);
indiv_optimal_k(rmv)=[];

for e=1:size(All_FC,4)
    for s=1:size(All_FC,3)
        sub=string(sub_list{s});
        if e==size(All_FC,4)
            if indiv_optimal_k(s)>0
                ts=dlmread(sprintf('%s/roi_ts_Schaef300_FIX_2P_%d_dicer.csv',sub,indiv_optimal_k(s)));
            else
                ts=dlmread(sprintf('%s/roi_ts_Schaef300_FIX_2P.csv',sub));
            end    
        else
            ts=dlmread(sprintf('%s/roi_ts_Schaef300_%s.csv',sub,string(process_names(e))));
        end
        [m,~]=size(ts);
        corr_mat = (ts * ts')/m;
        [evec,evals]=eig(corr_mat);
        evals=diag(evals);
        evals=sort(evals, 'descend');
        total_var=sum(sum(evals));
        var_explained_comp1(e,s) = evals(1,:)/total_var;
        disp(sub)
    end
end

%% Save relevant Results.
for i=1:length(preprocs)
    QC_FC(i,:)=All_QC(i).QCFC;
    dist_dep(i)=All_QC(i).DistDep;
end

save('Results/QC/QC_FC_corrs.mat', 'QC_FC');
save('Results/QC/distance_dependence.mat', 'dist_dep');
dlmwrite('Results/VE1/VE1_PC1.csv', var_explained_comp1)

