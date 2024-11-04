function check_erosions(anat_loc, orig_csf)

cd(anat_loc)

gunzip('wm_e5.nii.gz')
wm=niftiread('wm_e5.nii');
hdr=niftiinfo('wm_e5.nii');
i=5;
if sum(wm(:)>0)<5
    disp('WM mask has too few voxels')
    breakloop=0;
    while breakloop==0
        i=i-1;
        wm_name=sprintf('wm_e%u.nii',i);
        gunzip(wm_name+".gz")
        wm=niftiread(wm_name);
        hdr=niftiinfo(wm_name);
        if sum(wm(:)>0)>=5
            fprintf(1,'New WM Erosion is %u \n',i)
            breakloop=1;
        end
    end
end
wmmask=wm;
niftiwrite(wmmask,'final_wmmask.nii',hdr);
gzip('final_wmmask.nii')

csf=niftiread('csf_e2.nii.gz');
hdr=niftiinfo('csf_e2.nii.gz');
if sum(csf(:)>0)<5
    disp('mask has too few voxels, checking previous erosion')
    csf=niftiread('csf_e1.nii.gz');
    hdr=niftiinfo('csf_e1.nii.gz');
    if sum(csf(:)>0)<5
        disp('CSF erosions had too few voxels, taking original CSF mask')
        csf=niftiread(orig_csf);
        hdr=niftiinfo(orig_csf);
    elseif sum(csf(:)>0)>=5
        disp('First CSF erosion has enough voxels')
    end
end
csfmask=csf;
niftiwrite(csfmask,'final_csfmask.nii',hdr);
gzip('final_csfmask.nii')

end

