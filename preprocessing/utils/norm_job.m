function norm_job(rsdir, subname, dynnum, deform)

matlabbatch{1}.spm.spatial.normalise.write.subj.def = {deform};
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

matlabbatch{1}.spm.spatial.normalise.write.subj.resample=data1;
%%
matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = [-90 -126 -72
                                                          90 90 108];
matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = [2 2 2];
matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 4;
matlabbatch{1}.spm.spatial.normalise.write.woptions.prefix = 'w';

spm_jobman('initcfg')
spm_jobman('run',matlabbatch)
return
