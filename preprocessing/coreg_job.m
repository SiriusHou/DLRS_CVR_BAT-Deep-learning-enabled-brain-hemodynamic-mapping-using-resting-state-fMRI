function coreg_job(target,source,other)
% coregister source to target and apply transform matrix to source and 
% varargin
% dependent package - spm12

% test
% other = {'Y:\data\vcid\CVR_1001\meanhlu_vcid1001_cvr-d0281.img';'Y:\data\vcid\CVR_1001\rhlu_vcid1001_cvr-d0281.img'};

all_list   = {[source ',1']};

if nargin > 2
    num_otr     = length(other);
    for ii = 1:num_otr
        p_tmp   = spm_vol(other{ii});
        num_tmp = length(p_tmp);
        all_list = [all_list;strcat(repmat({[other{ii},',']},num_tmp,1),strsplit(num2str(1:num_tmp))')]; % n*1 cell
    end
end


matlabbatch{1}.spm.spatial.coreg.estimate.ref = {target};
matlabbatch{1}.spm.spatial.coreg.estimate.source = {[source ',1']};
matlabbatch{1}.spm.spatial.coreg.estimate.other = all_list; 
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];

spm('defaults','fmri'); 
spm_jobman('initcfg')
spm_jobman('run',matlabbatch)
return


