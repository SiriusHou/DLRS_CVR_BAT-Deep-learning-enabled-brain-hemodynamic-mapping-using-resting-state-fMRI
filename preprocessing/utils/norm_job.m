function norm_job(rsdir, subname, dynnum, tpm_file)

norm_est_file = spm_select('FPList', rsdir, ['^mean', subname, '.img']);
matlabbatch{1}.spm.spatial.normalise.estwrite.subj.vol = {norm_est_file};
for i=1:dynnum
    if i<10
        numstr=['00',num2str(i)];
    elseif i<100
        numstr=['0',num2str(i)];
    else 
        numstr=num2str(i);
    end
    
    boldfile=spm_select('FPlist', rsdir, 'img', ['^ar', subname, '-', numstr, '-001']);
    %disp(boldfile)
    data1{i,1}=[deblank(boldfile(1,:)),',1'];
    
end

matlabbatch{1}.spm.spatial.normalise.estwrite.subj.resample=data1;
%%
matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.biasreg = 0.0001;
matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.biasfwhm = 60;
matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.tpm = fullfile(tpm_file, filesep, "TPM.nii");
matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.affreg = 'mni';
matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.reg = [0 0.001 0.5 0.05 0.2];
matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.fwhm = 0;
matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.samp = 3;
matlabbatch{1}.spm.spatial.normalise.estwrite.woptions.bb = [-90 -126 -72
                                                          90 90 108];
matlabbatch{1}.spm.spatial.normalise.estwrite.woptions.vox = [2 2 2];
matlabbatch{1}.spm.spatial.normalise.estwrite.woptions.interp = 4;
matlabbatch{1}.spm.spatial.normalise.estwrite.woptions.prefix = 'w';

spm_jobman('initcfg')
spm_jobman('run',matlabbatch)
return
