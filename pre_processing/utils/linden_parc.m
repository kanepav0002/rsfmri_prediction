function linden_parc(func_folder, FC_out_name)
    
    func_folder=string(func_folder);
    fsl_ts=dlmread(func_folder+'/fsl_parc_ts.txt');
    %fsl_ts=zscore(fsl_ts);
    gm=dlmread(func_folder+'/fsl_gm_prob.txt');
    
    for i=1:size(fsl_ts,2)
        roi_TS(:,i)= fsl_ts(:,i)/gm(i);
    end

    FC=corr(roi_TS);
    FC=atanh(FC);
    
    func_folder=char(func_folder);
    sub_folder=func_folder; %(1:end-5);
    FC_out=(sub_folder+"/FC_"+FC_out_name);
    dlmwrite(FC_out,FC);
    roi_ts_out=(sub_folder+"/roi_ts_"+FC_out_name);
    dlmwrite(roi_ts_out,fsl_ts);
    gm_ts_out=(sub_folder+"/gmw_roi_ts_"+FC_out_name);
    dlmwrite(gm_ts_out,roi_TS);
end