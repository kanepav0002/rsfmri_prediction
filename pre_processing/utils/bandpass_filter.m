function bandpass_filter(func_folder, file_name, brain_mask, TR, lp, hp, repo_directory)
% This function bandpass filters a 4D niftiimage, needs a brain mask 3D
% file so that only brain voxels are filtered.
rest_loc=sprintf('%s/utils/REST-master',repo_directory);
addpath(genpath(rest_loc))

TR=str2double(TR);
lp=str2double(lp);
hp=str2double(hp);

to_load=append(func_folder,'/',file_name);
data=niftiread(to_load);
info=niftiinfo(to_load);
mask=niftiread(brain_mask);

% Keep only brain voxels for filtering
mask=reshape(mask,[],1);
mask=logical(mask);

dim=size(data);
data=reshape(data,[],dim(4));
data=data';
numVoxels=size(data,2);
time=size(data,1);
data=data(:,mask);

[filtered_dat] = rest_IdealFilter(data,TR,[lp hp]);

% Put filtered brain back into full image and reshape back to 4D
data_out=zeros(time,numVoxels);
data_out(:,mask)=filtered_dat;
data_out=data_out';
data_out=reshape(data_out,dim);
data_out=single(data_out);

name_out=extractBefore(file_name,".");
name_out=append(name_out,'_bpass.nii');
file_out=append(func_folder,'/',name_out);
niftiwrite(data_out,file_out,info)
gzip(file_out);
delete(file_out);
end
