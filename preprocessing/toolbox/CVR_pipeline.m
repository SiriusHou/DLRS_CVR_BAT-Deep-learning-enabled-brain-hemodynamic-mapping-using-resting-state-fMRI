function outDir = CVR_pipeline(commandLocation, paras)
% CVR_pipeline.m
% clear; close all; fclose all; clc;
% restoredefaultpath;
% addpath(genpath(['D:\cvrPipelineConstructionSite\'...
%     'peiyingYangStitchedPipeline_v5_git\toolbox']));
% addpath(genpath(['D:\cvrPipelineConstructionSite\'...
%     'peiyingYangStitchedPipeline_v5_git\spm12']));
% addpath(genpath(['D:\cvrPipelineConstructionSite\'...
%     'peiyingYangStitchedPipeline_v5_git\jsonlab']));
% commandLocation = ['D:\cvrPipelineConstructionSite\'...
%     'peiyingYangStitchedPipeline_v5_git\'...
%     'commands'];
% paras =  ['D:\cvrPipelineConstructionSite\' ...
%     'peiyingYangStitchedPipeline_v5_git\CO2O2_data\cvr_parameters_noMPR.json'];
% paras =  ['D:\cvrPipelineConstructionSite\peiyingYangStitchedPipeline_v5_git\vcid_1003_testEtCO2Only\cvr_parameters.json'];
% paras =  ['D:\cvrPipelineConstructionSite\peiyingYangStitchedPipeline_v5_git\CO2O2_data\cvr_parameters_noMPR.json'];
% paras =  ['D:\cvrPipelineConstructionSite\peiyingYangStitchedPipeline_v5_git\CO2O2_data\cvr_parameters_noMprSlc.json'];



%% Inputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load all user-input and other parameters from .json file
cvr_paras = loadjson(paras);

% get BOLD scan file name (.hdr file)                       - user input
boldHdrFileName = cvr_paras.userInput.boldHdrFileName;
% get BOLD scan file name (.img file)                       - user input
boldImgFileName = cvr_paras.userInput.boldImgFileName;
% get CO2 recording file name (.txt file)                   - user input
co2TraceFileName = cvr_paras.userInput.co2TraceFileName;
% get output of T1 Multiatlas Segmentation (.zip file)      - user input
multiatlasFileName = cvr_paras.userInput.multiatlasFileName;
% get the sample rate of the co2 recording device (in Hz)   - user input
sample_rate_co2 = cvr_paras.userInput.sample_rate_co2;
% get time of repitition (in seconds)                       - user input
TR = cvr_paras.userInput.TR;
% get order of slices (.txt file)                           - user input
sliceOrderFile = cvr_paras.userInput.sliceOrderFile;

% get the interpolation rate of the co2 envelope (in Hz)    - private
envelope_interp_rate = cvr_paras.private.envelope_interp_rate;
% get FWHM of gaussian smoothing kernal (in mm)             - private
SmoothFWHMmm = cvr_paras.private.SmoothFWHMmm;
% get etco2 peak picking control value                      - private
select = cvr_paras.private.select;
% get brain mask file name                                  - private
brainMaskName = cvr_paras.private.brainMaskName;
% get fixed delay range for voxelwise and ROIwise shift     - private
fixedDelayRange(1) = cvr_paras.private.fixedDelayRangeMin;
fixedDelayRange(2) = cvr_paras.private.fixedDelayRangeMax;

% project main path                                         - directory
mainPath = cvr_paras.Dir.mainPath;
% workspace file name                                       - directory
workDirFileName = cvr_paras.Dir.workDirFileName;
% output file name                                          - directory
outputFileName = cvr_paras.Dir.outputFileName;



%% General Setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get directory of TPM.nii
imgTpmDir = [commandLocation filesep 'TPM.nii'];

% if there is no mprage data, set mprageFlag off to disable all mprage
% processing
mprageFlag = 1;
if isempty(multiatlasFileName)
    mprageFlag = 0;
end

% create workDir in mainPath, which will contain all interim files
workDir = [mainPath filesep workDirFileName];
if exist(workDir,'dir')
    rmdir(workDir,'s');
end
mkdir(workDir)

% create outDir in mainPath, which will contain the outputs of the pipeline
outDir = [mainPath filesep outputFileName];
if exist(outDir,'dir')
    rmdir(outDir,'s');
end
mkdir(outDir)



% %% General Setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % load all user-input and other parameters from .json file
% cvr_paras = loadjson(paras);
% 
% % project main path
% mainPath = cvr_paras.Dir.mainPath;
% 
% % create workDir in mainPath, which will contain all interim files
% workDir = [mainPath filesep 'workspace'];
% if exist(workDir,'dir')
%     rmdir(workDir,'s');
% end
% mkdir(workDir)
% 
% % create outDir in mainPath, which will contain the outputs of the pipeline
% outDir = [mainPath filesep 'output'];
% if exist(outDir,'dir')
%     rmdir(outDir,'s');
% end
% mkdir(outDir)
% 
% 
% 
% %% Inputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % get BOLD scan file name (.hdr file)                       - user input
% boldHdrFileName = cvr_paras.userInput.boldHdrFileName;
% % get BOLD scan file name (.img file)                       - user input
% boldImgFileName = cvr_paras.userInput.boldImgFileName;
% % get CO2 recording file name (.txt file)                   - user input
% co2TraceFileName = cvr_paras.userInput.co2TraceFileName;
% % get output of T1 Multiatlas Segmentation (.zip file)      - user input
% multiatlasFileName = cvr_paras.userInput.multiatlasFileName;
% % get the sample rate of the co2 recording device (in Hz)   - user input
% sample_rate_co2 = cvr_paras.userInput.sample_rate_co2;
% % get time of repitition (in seconds)                       - user input
% TR = cvr_paras.userInput.TR;
% % get order of slices (.txt file)                           - user input
% sliceOrderFile = cvr_paras.userInput.sliceOrderFile;
% 
% % get the interpolation rate of the co2 envelope (in Hz)    - private
% envelope_interp_rate = cvr_paras.private.envelope_interp_rate;
% % get FWHM of gaussian smoothing kernal (in mm)             - private
% SmoothFWHMmm = cvr_paras.private.SmoothFWHMmm;
% % get etco2 peak picking control value                      - private
% select = cvr_paras.private.select;
% % get brain mask file name                                  - private
% brainMaskName = cvr_paras.private.brainMaskName;
% % get fixed delay range for voxelwise and ROIwise shift     - private
% fixedDelayRange(1) = cvr_paras.private.fixedDelayRangeMin;
% fixedDelayRange(2) = cvr_paras.private.fixedDelayRangeMax;
% 
% % get directory of TPM.nii                                  - directory
% imgTpmDir = [commandLocation filesep 'TPM.nii'];
% 
% % if there is no mprage data, set mprageFlag off to disable all mprage
% % processing
% mprageFlag = 1;
% if isempty(multiatlasFileName)
%     mprageFlag = 0;
% end
% 
% 

%% Relocate Input Files to workDir %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

copyfile([mainPath filesep boldHdrFileName],workDir,'f')
copyfile([mainPath filesep boldImgFileName],workDir,'f')
copyfile([mainPath filesep co2TraceFileName],workDir,'f')
if mprageFlag
copyfile([mainPath filesep multiatlasFileName],workDir,'f')
end
% sliceOrderFile is an optional input
if ~isempty(sliceOrderFile)
    copyfile([mainPath filesep sliceOrderFile],workDir,'f')
end



%% Process Additional Values %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get name of bold scan file without .img or .hdr extension
boldname = boldImgFileName(1:end-4);

% get bold scan's matrix size, matrix center, and voxel resolution (in mm)
fn_swrbold  = [workDir filesep boldImgFileName];
[boldvols, matsize, matresol] = read_hdrimg(fn_swrbold);
matcenter = matsize./2;

% get duration of bold scan
dynnum = size(boldvols,4);
bold_scandur = TR * dynnum;

if mprageFlag

% unzip segemented t1 info and return it's path infomation
[t1Path] = zach_unzipT1(workDir,multiatlasFileName);
% % get name of T1-weighted image file
% name_T1 = multiatlasFileName(1:end-4);
% 
% % create path to unzipped data
% t1Path = [workDir filesep name_T1];
% if exist(t1Path,'dir')
%     rmdir(t1Path,'s')
% end
% mkdir(t1Path)
% 
% % unzip anatomical T1-weighted image data (dump it into workDir)
% unzip([workDir filesep name_T1 '.zip'],t1Path)
end

% Get .mat info from first dynamic (to be used during coregistration)
P = spm_vol(fn_swrbold);
matInfo = P(1).mat;



%% Global Parameter Setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd(workDir);
spm_get_defaults;
global defaults;
defaults.mask.thresh = 0;



%% Slice Timing Correction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% skip this step if sliceOrderFile is omitted
if ~isempty(sliceOrderFile)

    % get original bold scan
    P = spm_select('List', workDir, ['^' boldname '.img']);

    % set TA and refSlice
    scanInfo = spm_vol(P);
    nslices = scanInfo(1).dim(3); % use first dynamic to get number of slices
    TA = TR-(TR/nslices);
    refSlice = 1;

    % get order of slices
    fid = fopen([workDir filesep sliceOrderFile]);
    sliceOrder = textscan(fid,'%f');
    sliceOrder = sliceOrder{1};
    fclose(fid);

    % get timing for correction
    timing(2)=TR-TA;
    timing(1)=TA/(nslices-1);

    % correct slice timing (produces a new bold scan with the filename prefixed
    % with 'a')
    spm_slice_timing(P, sliceOrder, refSlice, timing);

    % rename files so boldname reflects slice-time-corrected scan
    movefile([workDir filesep boldname '.img'],...
        [workDir filesep 'b4SliceTimeCorrection_' boldname '.img'])
    movefile([workDir filesep boldname '.hdr'],...
        [workDir filesep 'b4SliceTimeCorrection_' boldname '.hdr'])
    movefile([workDir filesep 'a' boldname '.img'],...
        [workDir filesep boldname '.img'])
    movefile([workDir filesep 'a' boldname '.hdr'],...
        [workDir filesep boldname '.hdr'])
end



%% Realign BOLD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get original bold scan
P = cell(1,1);    
P{1}   = spm_select('List',workDir,['^' boldname '.*img']);

% get the bold scan's data
V       = spm_vol(P); 
V       = cat(1,V{:});

% realign bold scan
disp(sprintf(['realigning ' boldname]));
FlagsC = struct('quality',defaults.realign.estimate.quality,...
    'fwhm',5,'rtm',0);
spm_realign(V, FlagsC);

% reslice bold scan
which_writerealign = 2;
mean_writerealign = 1;
FlagsR = struct('interp',defaults.realign.write.interp,...
    'wrap',defaults.realign.write.wrap,...
    'mask',defaults.realign.write.mask,...
    'which',which_writerealign,'mean',mean_writerealign); 
spm_reslice(P,FlagsR);



%% Smooth BOLD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get realigned and resliced bold scan (prefixed with 'r')
P = cell(1,1);
P{1}   = spm_select('List',workDir,['^r' boldname '.*img']);

% get the scan's data
V       = spm_vol(P);
V       = cat(1,V{:});

% smooth scan (creates a 4D smoothed scan prefixed with 's' and the FWHM of
% the gaussian kernel) 
disp(sprintf(['smoothing ' boldname]));
[pth,nam,ext] = fileparts(V(1).fname);
fnameIn       = fullfile(pth,[nam ext]);
fname         = fullfile(pth,['s' int2str(SmoothFWHMmm) nam ext]);          
spm_smooth(fnameIn,fname,SmoothFWHMmm);

% clear V for further use
clear V;



%% Find EtCo2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get filename for the upper envelope of the CO2 trace
envelopeCo2FilePath = zach_standaloneFunc_getEtco2Envelope_v2(workDir,...
    co2TraceFileName,sample_rate_co2,envelope_interp_rate,select);



%% Shift EtCO2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get average smoothed bold scan (prefixed with 'mean')
meanimg = spm_select('List',workDir,['^mean' boldname '.*img']);

% get brain masks using average bold scan
[brnmsk, brnmsk_clcu] = bold_getBrainMask(imgTpmDir,meanimg);

% save masks as 'BrainMask.img' and 'BrainMask_calulation.img'
write_ANALYZE(brnmsk,[brainMaskName '.img'],matsize,matresol,1,4,0,...
    matcenter);
write_ANALYZE(brnmsk_clcu,[brainMaskName '_calulation.img'],matsize,...
    matresol,1,4,0,matcenter);
bold_maskPath = [workDir filesep brainMaskName '.img'];

% get realigned and resliced bold scan and its respective data
P = cell(1,1);
P{1} =  spm_select('List',workDir,['^r' boldname '.*img']);
V       = spm_vol(P);
V       = cat(1,V{:});

% find the average bold signal for the whole brain inside the
% calculation brian mask (this is a 1D signal)
avgBold = zeros(length(V),1);
for i=1:length(V)
    vols = spm_read_vols(V(i));
    vols(isnan(vols)) = 0;
    img = vols/V(i).pinfo(1);
    avgBold(i) = mean(img(brnmsk_clcu>0));
end

% save the whole-brain bold signal
name_avgboldPath=strcat(workDir,filesep,'AvgBOLD_WB_ns.txt');
save(name_avgboldPath, 'avgBold','-ascii');

% get the whole-brain bold signal
avgBold=textread(name_avgboldPath,'%f');

% get the co2 envelope signal (timestamps in 1st column, mmHgCo2 in 2nd)
r1 = textread(envelopeCo2FilePath,'%f','delimiter',',');
etco2timecourse = reshape(r1,2,length(r1)/2)';

% find the optimal delay between the whole-brain bold signal and the co2
% envelope (essentially done by finding the lowest residual values between
% the two shifted curves (after linear fitting))

% delay range in seconds (it is assumed that the EtCO2 curve will cover a
% much larger range than the BOLD; therefore, the shift of the EtCO2 curve
% will only be negative (to the left)) 
minRange = -abs((length(etco2timecourse)./envelope_interp_rate) -...
     (length(avgBold).*TR));
maxRange = 0;
delayrange = [minRange maxRange];

% finds the optimal delay in the delay range between the two curves
[optDelay, optEtCO2] = cvr_func_findCO2delay(avgBold, TR,...
    etco2timecourse, delayrange,1,1,1,outDir);

% repeat process once with +/- 10s delays with 0.1s iteration bewteen
% delays to approach (zoom in on) the optimal delay
delayrange = [optDelay-5, optDelay+5];
[optDelay, optEtCO2] = cvr_func_findCO2delay(avgBold, TR,...
    etco2timecourse, delayrange,1,0,0.1,outDir);

co2delay=optDelay;

% save delays
filename = fullfile(workDir, 'Sync_EtCO2_timecourse.txt');
save(filename,'optEtCO2','-ascii'); 
filename = fullfile(workDir, 'EtCO2_BOLD_delay.txt');
save(filename,'co2delay','-ascii');

% use delay to generate an EtCO2 file with modified timestamps,
% where t=0 is at the beginning of the EtCO2 synced with the global BOLD
% signal
boldSyncedEtco2Path = zach_syncEtco2WithBold(envelopeCo2FilePath,...
    co2delay);

% generate figure comparing whole brain BOLD signal with boldSyncedEtco2
zach_co2BoldfigureCompGenerator(boldSyncedEtco2Path,name_avgboldPath,...
    outDir,TR);


%% Global-shift CVR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% find top 25% and get delta (top-bttm) %%%%%%%%%%%%%%%%
% use mean/min etco2 in future %%%%%%%%%%%%%%%%%%%%%%%%

% find and save CVR map using global-shift method
cvrdir_g = [workDir filesep 'CVR_globalshift'];mkdir(cvrdir_g);
brnmsk=loadimage([workDir filesep brainMaskName '.img'],1);
[EtCO2_mean,EtCO2_min,EtCO2_max] = CVR_mapping_spm_GLM(workDir,boldname,...
    TR,SmoothFWHMmm,brnmsk,0,cvrdir_g);

% % find and save relative CVR map as well (prefixes files with 'rel_')
% globalShiftCvrFileName = ['HC_CVRmap_s' int2str(SmoothFWHMmm)];
% relMap = zachCVR_generateRelativeCvrMap(cvrdir_g,globalShiftCvrFileName,...
%     brnmsk);
% fn_relCvr = [cvrdir_g filesep 'rel_' globalShiftCvrFileName '.img'];
% write_hdrimg(relMap,fn_relCvr,matresol,16,1);

% convert nii to img/hdr files
niis = {[cvrdir_g filesep 'beta_0001.nii'];
        [cvrdir_g filesep 'beta_0002.nii'];
        [cvrdir_g filesep 'beta_0003.nii']};
for inii = 1:3
    p_nii = spm_vol(niis{inii});
    v_nii = spm_read_vols(p_nii);
    p_img = p_nii;
    p_img.fname = strrep(niis{inii},'.nii','.img');
    spm_write_vol(p_img,v_nii);
end

% get paths to globalshift beta values
bold_gBeta1Path = [cvrdir_g filesep 'beta_0001.img'];
bold_gBeta2Path = [cvrdir_g filesep 'beta_0002.img'];
bold_gBeta3Path = [cvrdir_g filesep 'beta_0003.img'];



%% Voxel-shift Beta Maps %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% find beta maps, EtCO2_mean, and EtCO2_min (used to find CVR map)and BAT
% map using the voxel-shift method
cvrdir_v = [workDir filesep 'CVR_voxelshift_etco2'];mkdir(cvrdir_v);
voxelwiseResult = CVR_mapping_voxelwise_GLM(workDir, boldname, TR,...
    SmoothFWHMmm, fixedDelayRange, brnmsk, boldSyncedEtco2Path, cvrdir_v);

% get paths to BAT map and voxelshift beta values
bold_batPath = voxelwiseResult{1};
bold_vBeta1Path = voxelwiseResult{2};
bold_vBeta2Path = voxelwiseResult{3};
bold_vBeta3Path = voxelwiseResult{4};

% voxelShiftCvrFileName = [boldname '_s' int2str(SmoothFWHMmm) '_CVR'];
% fn_cvr = [cvrdir_v filesep voxelShiftCvrFileName '.img'];
% beta1 = spm_select('FPList',cvrdir_v,'_beta_0001.img$');
% beta3 = spm_select('FPList',cvrdir_v,'_beta_0003.img$');



%% Skip Further Processing if MPRAGE Data is Omitted %%%%%%%%%%%%%%%%%%%%%%
if mprageFlag



%% Interum Value Setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get mean bold image file name
meanBoldName = ['mean',boldImgFileName(1:end-4)];

% get T1 file name
t1ImgSizeFile = spm_select('List',t1Path,'.*[^mni].imgsize$');
t1ImgName = t1ImgSizeFile(1:end-8);

% get skull-stripped MPRAGE brain scan
mpr_brain = ASL_mprageSkullstrip(t1Path, t1ImgName);

% grab smoothed BOLd for interum value return
smoothBoldPath = spm_select('FPList',workDir,...
    ['^s' int2str(SmoothFWHMmm) 'r' '.*\.img']); 



%% Coregistration %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Split 4D boldScan .img and .hdr into 3D .img and .hdr files; put files
% into directory named after boldImgCoregName
boldImgPath = [workDir filesep 'r' boldImgFileName];
[~,boldImgName] = fileparts(boldImgPath);
dumpDir = [workDir filesep 'r'];
if exist(dumpDir,'dir')
    rmdir(dumpDir,'s')
end
mkdir(dumpDir); addpath(dumpDir);
boldScans3d = zach_split4DScan(boldImgPath,dumpDir);
boldScans3dDir = dumpDir;

% do the same for the smoothed BOLD
[~,boldImgName] = fileparts(smoothBoldPath);
dumpDirS = [workDir filesep 's'];
if exist(dumpDirS,'dir')
    rmdir(dumpDirS,'s')
end
mkdir(dumpDirS); addpath(dumpDirS);
boldScansS3d = zach_split4DScan(smoothBoldPath,dumpDirS);
boldScansS3dDir = dumpDirS;

% I don't know what this does but I'm scared of removing it
% bold_maskPath = [workDir filesep brainMaskName '.img'];
brnmsk = spm_read_vols(spm_vol(bold_maskPath));
brnmsk = brnmsk > 0.5;

% reset the .mat of each hdr/img so that they are aligned
fl1 = {   [workDir filesep meanBoldName '.img'];
          bold_gBeta1Path;
          bold_gBeta2Path;
          bold_gBeta3Path;
          bold_vBeta1Path;
          bold_vBeta2Path;
          bold_vBeta3Path;
          bold_batPath;
          bold_maskPath};
for jj = 1:length(fl1)
    p_img = spm_vol(fl1{jj});
    v_img = spm_read_vols(p_img);
    p_img.mat = matInfo; % acquired at Process Additional Values subheading
    spm_write_vol(p_img,v_img);
end

% coregister maps to skull-stripped MPRAGE
target = mpr_brain;
source = [workDir filesep meanBoldName '.img'];
other = { bold_gBeta1Path;
          bold_gBeta2Path;
          bold_gBeta3Path;
          bold_vBeta1Path;
          bold_vBeta2Path;
          bold_vBeta3Path;
          bold_batPath;
          bold_maskPath};
other = [other;boldScans3d];
other = [other;boldScansS3d];
asl_coreg12(target,source,other);

% get coregistered beta1 map paths
[pth,nme,ext] = fileparts(bold_gBeta1Path);
mpr_gBeta1Path = [pth filesep 'r' nme ext];
[pth,nme,ext] = fileparts(bold_vBeta1Path);
mpr_vBeta1Path = [pth filesep 'r' nme ext];

% get coregistered beta2 map paths
[pth,nme,ext] = fileparts(bold_gBeta2Path);
mpr_gBeta2Path = [pth filesep 'r' nme ext];
[pth,nme,ext] = fileparts(bold_vBeta2Path);
mpr_vBeta2Path = [pth filesep 'r' nme ext];

% get coregistered beta3 map paths
[pth,nme,ext] = fileparts(bold_gBeta3Path);
mpr_gBeta3Path = [pth filesep 'r' nme ext];
[pth,nme,ext] = fileparts(bold_vBeta3Path);
mpr_vBeta3Path = [pth filesep 'r' nme ext];

% get coregistered BAT map path
[pth,nme,ext] = fileparts(bold_batPath);
mpr_batPath = [pth filesep 'r' nme ext];

% get coregistered brain mask path
[pth,nme,ext] = fileparts(bold_maskPath);
mpr_maskPath = [pth filesep 'r' nme ext];

% move coregistered BOLD scan paths to new directory
mpr_3dBoldScansDir = [boldScans3dDir filesep 'rr'];
if exist(mpr_3dBoldScansDir,'dir')
    rmdir(mpr_3dBoldScansDir,'s')
end
mkdir(mpr_3dBoldScansDir); addpath(mpr_3dBoldScansDir);
filesToMove = spm_select('FPList',boldScans3dDir,'^rr');
for i = 1:size(filesToMove,1)
    path = filesToMove(i,:);
    movefile(path,mpr_3dBoldScansDir,'f')
end

% move coregistered smoothed BOLD scan paths to new directory
mpr_3dBoldScansSDir = [boldScansS3dDir filesep 'rr'];
if exist(mpr_3dBoldScansSDir,'dir')
    rmdir(mpr_3dBoldScansSDir,'s')
end
mkdir(mpr_3dBoldScansSDir); addpath(mpr_3dBoldScansSDir);
filesToMove = spm_select('FPList',boldScansS3dDir,'^rs');
for i = 1:size(filesToMove,1)
    path = filesToMove(i,:);
    movefile(path,mpr_3dBoldScansSDir,'f')
end



%% Convert to MNI Space %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% normalize MPR beta maps, MPR bat map, and MPR brain mask to MNI space
mni_gBeta1Path = zach_mniNormalize(mpr_gBeta1Path, commandLocation, t1Path);
mni_gBeta2Path = zach_mniNormalize(mpr_gBeta2Path, commandLocation, t1Path);
mni_gBeta3Path = zach_mniNormalize(mpr_gBeta3Path, commandLocation, t1Path);
mni_vBeta1Path = zach_mniNormalize(mpr_vBeta1Path, commandLocation, t1Path);
mni_vBeta2Path = zach_mniNormalize(mpr_vBeta2Path, commandLocation, t1Path);
mni_vBeta3Path = zach_mniNormalize(mpr_vBeta3Path, commandLocation, t1Path);
mni_batPath = zach_mniNormalize(mpr_batPath, commandLocation, t1Path);
mni_maskPath = zach_mniNormalize(mpr_maskPath, commandLocation, t1Path);

% Normalize coregistered bold scan files to MNI space
mni_boldDir = zach_mniNormalize3DFileDir(mpr_3dBoldScansDir,...
    commandLocation, t1Path);

% Normalize coregistered smoothed bold scan files to MNI space
mni_boldSDir = zach_mniNormalize3DFileDir(mpr_3dBoldScansSDir,...
    commandLocation, t1Path);

% % put coregistered beta maps, bat map, and brain mask into MNI space
% mniBeta1 = zach_mniNormalize(mpr_vBeta1, commandLocation, t1Path);
% mniBeta3 = zach_mniNormalize(mpr_vBeta3, commandLocation, t1Path);
% mniBat2mmPath = zach_mniNormalize(mpr_batPath, commandLocation,...
%     t1Path);



% %get MNI brainmask from first MNI bold scan
% boldmni1 = spm_select('FPList',mni_boldDir,'ts001_MNI.img$');
% [~,~,brnmskMNI] = asl_getBrainMask(imgTpmDir,boldmni1);
% p_out = spm_vol(boldmni1);
% p_out.fname = [workDir filesep 'BrainMask_mni.img'];
% mniMaskFile2mmPath = p_out.fname;
% p_out.dt = [4 0];
% spm_write_vol(p_out,brnmskMNI);

% end of MPRAGE-dependent processing
end


%% Calculate CVR Maps

% get bold space globalshift CVR maps                             TO OUTPUT
result1 = zach_CVR_calc_betamaps(...
    [outDir filesep 'BOLDspace_globalshift_CVR.img'],bold_gBeta1Path,...
    bold_gBeta3Path, EtCO2_mean, EtCO2_min, bold_maskPath);
bold_globalshiftCvrPath = result1{1};
rel_bold_globalshiftCvrPath = result1{2};

% get bold space voxelshift CVR maps                              TO OUTPUT
result2 = zach_CVR_calc_betamaps(...
    [outDir filesep 'BOLDspace_voxelshift_CVR.img'],bold_vBeta1Path,...
    bold_vBeta3Path, EtCO2_mean, EtCO2_min, bold_maskPath);
bold_voxelshiftCvrPath = result2{1};
rel_bold_voxelshiftCvrPath = result2{2};

if mprageFlag
% get mpr space gloablshift CVR maps                              TO OUTPUT
result3 = zach_CVR_calc_betamaps(...
    [outDir filesep 'MPRspace_globalshift_CVR.img'],mpr_gBeta1Path,...
    mpr_gBeta3Path, EtCO2_mean, EtCO2_min, mpr_maskPath);
mpr_globalshiftCvrPath = result3{1};
rel_mpr_globalshiftCvrPath = result3{2};

% get mpr space voxelshift CVR maps                               TO OUTPUT
result4 = zach_CVR_calc_betamaps(...
    [outDir filesep 'MPRspace_voxelshift_CVR.img'],mpr_vBeta1Path,...
    mpr_vBeta3Path, EtCO2_mean, EtCO2_min, mpr_maskPath);
mpr_voxelshiftCvrPath = result4{1};
rel_mpr_voxelshiftCvrPath = result4{2};

% get mni space globalshift CVR maps                              TO OUTPUT
result5 = zach_CVR_calc_betamaps(...
    [outDir filesep 'MNIspace_globalshift_CVR.img'],mni_gBeta1Path,...
    mni_gBeta3Path, EtCO2_mean, EtCO2_min, mni_maskPath);
mni_globalshiftCvrPath = result5{1};
rel_mni_globalshiftCvrPath = result5{2};
[pth,nme,ext] = fileparts(mni_globalshiftCvrPath);
ASL_downSampleMNI(mni_globalshiftCvrPath,[pth filesep nme '_2mm' ext])
[pth,nme,ext] = fileparts(rel_mni_globalshiftCvrPath);
ASL_downSampleMNI(rel_mni_globalshiftCvrPath,[pth filesep nme '_2mm' ext])

% get mni space voxelshift CVR maps                               TO OUTPUT
result6 = zach_CVR_calc_betamaps(...
    [outDir filesep 'MNIspace_voxelshift_CVR.img'],mni_vBeta1Path,...
    mni_vBeta3Path, EtCO2_mean, EtCO2_min, mni_maskPath);
mni_voxelshiftCvrPath = result6{1};
rel_mni_voxelshiftCvrPath = result6{2};
[pth,nme,ext] = fileparts(mni_voxelshiftCvrPath);
ASL_downSampleMNI(mni_voxelshiftCvrPath,[pth filesep nme '_2mm' ext])
[pth,nme,ext] = fileparts(rel_mni_voxelshiftCvrPath);
ASL_downSampleMNI(rel_mni_voxelshiftCvrPath,[pth filesep nme '_2mm' ext])
end

% % caluculate CVR map from beta maps
% cvrMap = zach_CVR_calc_betamaps(spm_read_vols(spm_vol(mpr_vBeta1)),...
%     spm_read_vols(spm_vol(mpr_vBeta3)), EtCO2_mean, EtCO2_min,...
%     mpr_maskPath);

% % write CVR map to hdr/img files
% p_out1  = spm_vol(mpr_vBeta1);
% p_out1.fname = [cvrdir_v filesep 'rcvrmap_voxshft.img'];
% voxelShiftCvrFilePathWithBetaCoreg = p_out1.fname;
% spm_write_vol(p_out1,cvrMap);

% % find and save the relative CVR map as well (prefixes files with 'rel_')
% [pth,nme] = fileparts(voxelShiftCvrFilePathWithBetaCoreg);
% relMap = zachCVR_generateRelativeCvrMap(pth,nme,...
%     spm_read_vols(spm_vol(mpr_maskPath)));
% fn_relCvr = [cvrdir_v filesep 'rel_' nme '.img'];
% write_hdrimg(relMap,fn_relCvr,matresol,16,1);

% % use beta maps in MNI space to find CVR map in MNI space
% cvrMapMNI = zach_CVR_calc_betamaps(spm_read_vols(spm_vol(mniBeta1)),...
%     spm_read_vols(spm_vol(mniBeta3)), EtCO2_mean, EtCO2_min,...
%     mniMaskFile2mmPath);
% p_out1  = spm_vol(mniBeta1);
% p_out1.fname = [cvrdir_v filesep 'rcvrmap_mni_2mm_voxshift.img'];
% cvrPathMni2mm = p_out1.fname;
% spm_write_vol(p_out1,cvrMapMNI);

% % find and save the relative CVR map as well (prefixes files with 'rel_')
% [pth,nme] = fileparts(cvrPathMni2mm);
% relMap = zachCVR_generateRelativeCvrMap(pth,nme,...
%     spm_read_vols(spm_vol(mpr_maskPath)));
% fn_relCvrMni = [cvrdir_v filesep 'rel_' nme '.img'];
% write_hdrimg(relMap,fn_relCvrMni,matresol,16,1);



%% Find ROI Values %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% find CVR values by ROI and send values      

% for partial scans, only allow values to be exported if 50% of voxels are
% present

% provide 289 rois as well

% proivide voxel count for mni-mprage rois and mni-masked mni-bold rois
% (this is how you will get percent coverage of present voxels)

% TO OUTPUT
if mprageFlag
zach_cvrMappingROIwise(mni_boldDir, mni_maskPath,...
    boldSyncedEtco2Path, TR, fixedDelayRange, t1Path, outDir)
end
% rename BAT to relative BAT



%% Send to Output %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% send globalshift BOLD space cvr map                             TO OUTPUT

% send globalshift MPRAGE space cvr map                           TO OUTPUT

% send globalshift MNI space cvr map                              TO OUTPUT

% send voxelshift BOLD space cvr map                              TO OUTPUT



% % send voxelshift MPRAGE space cvr map                            TO OUTPUT
% copyfile(voxelShiftCvrFilePathWithBetaCoreg,outDir)
% [pth,nme] = fileparts(voxelShiftCvrFilePathWithBetaCoreg);
% copyfile([pth filesep nme '.hdr'],outDir)
%     % relative
% copyfile(fn_relCvr,outDir)
% [pth,nme] = fileparts(fn_relCvr);
% copyfile([pth filesep nme '.hdr'],outDir)
% 
% % send voxelshift MNI space cvr map                               TO OUTPUT
% [pth,nme,ext] = fileparts(cvrPathMni2mm);
% copyfile([pth filesep nme ext],outDir)
% copyfile([pth filesep nme '.hdr'],outDir)
%     % relative
% [pth,nme,ext] = fileparts(fn_relCvrMni);
% copyfile([pth filesep nme ext],outDir)
% copyfile([pth filesep nme '.hdr'],outDir)


% send BOLD space BAT map                                         TO OUTPUT
copyfile(bold_batPath,[outDir filesep 'rel_BOLDspace_BAT.img'])
[pth,nme] = fileparts(bold_batPath);
copyfile([pth filesep nme '.hdr'],[outDir filesep 'rel_BOLDspace_BAT.hdr'])

if mprageFlag
% send MPRAGE space BAT map                                       TO OUTPUT
copyfile(mpr_batPath,[outDir filesep 'rel_MPRspace_BAT.img'])
[pth,nme] = fileparts(mpr_batPath);
copyfile([pth filesep nme '.hdr'],[outDir filesep 'rel_MPRspace_BAT.hdr'])

% send MNI space BAT map                                          TO OUTPUT
copyfile(mni_batPath,[outDir filesep 'rel_MNIspace_BAT.img'])
[pth,nme] = fileparts(mni_batPath);
copyfile([pth filesep nme '.hdr'],[outDir filesep 'rel_MNIspace_BAT.hdr'])
% (include 2x2x2 MNIspace BAT map)
mniBATname = [outDir filesep 'rel_MNIspace_BAT.img'];
[pth,nme,ext] = fileparts(mniBATname);
ASL_downSampleMNI(mniBATname,[pth filesep nme '_2mm' ext])
end

% [pth,nme,ext] = fileparts(mniBat2mmPath);
% copyfile([pth filesep nme ext],outDir)
% copyfile([pth filesep nme '.hdr'],outDir)



% send EtCO2 mean/max/min/globalBAT .txt                          TO OUTPUT
fid = fopen([outDir filesep 'EtCO2_Stats.txt'],'w');
fprintf(fid, ['Mean_EtCO2(mmHg)\tMaximal_EtCO2(top_25%%_avg)(mmHg)\t'...
    'Minimal_EtCO2(bttm_25%%_avg)(mmHg)\tGlobal_EtCO2_Shift(s)\n']);
fprintf(fid, [num2str(EtCO2_mean) '\t' num2str(EtCO2_max) '\t'...
    num2str(EtCO2_min) '\t' num2str(co2delay)]);
fclose(fid);


%% Identify Outputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

outputFileNames = dir(outDir);
outputFileNames = {outputFileNames.name};
for i = 1:length(outputFileNames)
    if strcmp(outputFileNames{i},'.') || strcmp(outputFileNames{i},'..')
        continue
    end
    movefile([outDir filesep outputFileNames{i}],...
        [outDir filesep boldname '_' outputFileNames{i}]);
end



%% Finish and cleanup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
disp('Processing Completed')



