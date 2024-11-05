function segment_anat(T1, repo_directory)
    
    T1=char(T1);
    linden_scripts=sprintf('%s/utils/linden_rsfmri',repo_directory);
    spm_dir=sprintf('%s/utils/spm8/matlab2015b.r6685/',repo_directory);
    addpath(genpath(linden_scripts))
    addpath(genpath(spm_dir))
    gunzip(T1);
    T1=T1(1:end-3);
    SegmentT1(T1, spm_dir);
end
