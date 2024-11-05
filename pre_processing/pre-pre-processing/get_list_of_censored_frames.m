%% Create Censored Frames list

% Inputs
data_dir="put the directory of the BIDS formatted data directory here";
sub_list=readcell('a csv/txt file that contains a 1xN list of subject IDs');
fd_file_name="the name of the FD file for each subject";

% script
sub_idx=1;
for s=1:length(sub_list)
    sub=string(sub_list{s});
    
    fd=dlmread(sprintf('%s/%s/func/%s',data_dir,sub,fd_file_name));  
    tmp_censor=fd<0.3;
    

    if sum(tmp_censor(1:5)) < 5
        tmp_censor(1:5)=0;
    end
    
    for i = 6:(length(tmp_censor)-6)
        curr_val=tmp_censor(i);
        before_sum=sum(tmp_censor(i-4:i));
        after_sum=sum(tmp_censor(i:i+4));
        middle_sum=(sum(tmp_censor(i-2:i)))+(sum(tmp_censor(i+1:i+2)));
        
        if before_sum <5 && after_sum <5 && middle_sum < 5
            tmp_censor(i) = 0;
        end
    end
    
    n_frames=length(tmp_censor);
    if sum(tmp_censor(n_frames-5:end)) <5
        tmp_censor(n_frames-5:end)=0;
    end
    
    file_out=sprintf('%s/%s/func/censored_frames.txt',data_dir,sub);
    dlmwrite(file_out,tmp_censor)
    disp(sub)
end
    
