function smooth(func_file, func_folder, repo_directory)
    
    linden_scripts=sprintf('%s/utils/linden_rsfmri',repo_directory);
    spm_dir=sprintf('%s/utils/spm8/matlab2015b.r6685/',repo_directory);
    addpath(genpath(linden_scripts))
    addpath(genpath(spm_dir))
    
    func=char(func_folder+"/"+func_file);
    gunzip(func)
    [hdr,data]=read(func(1:end-3));
    dim=size(data);
    gunzip(func);
    func=(func(1:end-3));
    cd(func_folder)
    func_file=char(func_file);
    
    SmoothEPI(func_file(1:end-3),6,dim(4));
    func_file=("s"+func_file(1:end-3));
    gzip(func_file) 
    
end
