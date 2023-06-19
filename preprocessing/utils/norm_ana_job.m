function norm_ana_job(target, deform)
matlabbatch{1}.spm.spatial.normalise.write.subj.def = {deform};

matlabbatch{1}.spm.spatial.normalise.write.subj.resample = {target};
matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = [-90 -126 -72
                                                          90 90 108];
matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = [2 2 2];
matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 4;
matlabbatch{1}.spm.spatial.normalise.write.woptions.prefix = 'w';
spm_jobman('initcfg')
spm_jobman('run',matlabbatch)
return