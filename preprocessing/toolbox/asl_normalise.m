function asl_normalise(source,varargin)
% normalize source to nmi tpm, and apply transform matrix to source and
% varargin
% dependent package - spm12
% input:
% 	imgtpm, full path of the template
% 	source, full path of the source image to be normalized
% 	varargin, 1*n cell matrix in which each cell element contains a string of full path

imgtpm = 'C:\Users\yli199\Documents\MATLAB\spm12\tpm\TPM.nii';

spm('defaults', 'FMRI');

paras_normalise.subj.vol = {[source ',1']};
if iscell(varargin{1})
    paras_normalise.subj.resample = varargin{1};
else
    paras_normalise.subj.resample = strcat(varargin', ',1');
end
paras_normalise.eoptions.biasreg = 0.0001;
paras_normalise.eoptions.biasfwhm = 60;
paras_normalise.eoptions.tpm = {imgtpm};
paras_normalise.eoptions.affreg = 'mni';
paras_normalise.eoptions.reg = [0 0.001 0.5 0.05 0.2];
paras_normalise.eoptions.fwhm = 0;
paras_normalise.eoptions.samp = 3;
paras_normalise.woptions.bb = [-90 -126 -72
                                90   90 108]; % 2*2*2 mni
paras_normalise.woptions.vox = [2 2 2];
% paras_normalise.woptions.interp = 1; % trilinear
paras_normalise.woptions.interp = 4; % default: 4th b-spline
paras_normalise.woptions.prefix = 'w';

matlabbatch = [];
matlabbatch{1}.spm.spatial.normalise.estwrite = paras_normalise;
spm_jobman('run', matlabbatch);

