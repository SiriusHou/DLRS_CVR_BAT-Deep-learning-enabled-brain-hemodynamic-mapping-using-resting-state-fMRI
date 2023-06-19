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

spm('defaults','fmri'); spm_jobman('initcfg');

coregparas = struct();

coregparas.ref = {target};
coregparas.source = {[source ',1']};
coregparas.other = all_list; 
coregparas.eoptions.cost_fun = 'nmi';
coregparas.eoptions.sep = [4 2];
coregparas.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
coregparas.eoptions.fwhm = [7 7];


matlabbatch = [];
matlabbatch{1}.spm.spatial.coreg.estwrite = coregparas;
spm_jobman('run',matlabbatch);
