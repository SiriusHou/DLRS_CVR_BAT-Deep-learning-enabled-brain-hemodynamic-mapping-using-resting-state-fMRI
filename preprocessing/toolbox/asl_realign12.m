function asl_realign12(imgarray,varargin)
% realign image arrays to do motion correction
% dependent package - spm12

% get image file names
V       = spm_vol(imgarray);
vnum    = length(V);
namearray = strcat(repmat({[imgarray ',']},vnum,1),strsplit( num2str(1:vnum))');

realignparas = struct();

spm('defaults','fmri'); spm_jobman('initcfg');

realignparas.data = {namearray}'; % realignparas.data = namearray;
realignparas.eoptions.quality = 0.9;
realignparas.eoptions.sep = 4;
realignparas.eoptions.fwhm = 5;
realignparas.eoptions.rtm = 1;
realignparas.eoptions.interp = 2;
realignparas.eoptions.wrap = [0 0 0];
realignparas.eoptions.weight = '';
realignparas.roptions.which = [2 1];
realignparas.roptions.interp = 1; % trilinear
realignparas.roptions.wrap = [0 0 0];
realignparas.roptions.mask = 1; % apply mask or problem when do std
realignparas.roptions.prefix = 'r';

matlabbatch = [];
matlabbatch{1}.spm.spatial.realign.estwrite = realignparas;

spm_jobman('run',matlabbatch);

