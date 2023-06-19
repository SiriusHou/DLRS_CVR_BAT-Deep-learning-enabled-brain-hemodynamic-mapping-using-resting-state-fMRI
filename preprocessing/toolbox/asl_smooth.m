function asl_smooth(imgarray,varargin)
% realign image arrays to do motion correction
% dependent package - spm12

% get image file names
V           = spm_vol(imgarray);
vnum        = length(V);
namearray   = strcat(repmat({[imgarray ',']},vnum,1),strsplit( num2str(1:vnum))');
sw          = ones(1,3) * varargin{1};

spm('defaults','FMRI'); spm_jobman('initcfg');

smoothparas = struct();

smoothparas.data = namearray;
smoothparas.fwhm = sw;
smoothparas.dtype = 0;
smoothparas.im = 0;
smoothparas.prefix = ['s' num2str(varargin{1})];

matlabbatch = [];
matlabbatch{1}.spm.spatial.smooth = smoothparas;

spm_jobman('run',matlabbatch);

