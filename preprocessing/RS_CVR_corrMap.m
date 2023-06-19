function RS_CVR_corrMap(wkdir, subname, code_directory)

warning off;
% Global Parameter Setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mni_resolution = [2, 2, 2];
mni_dimension = [91, 109, 91];
mni_type = 16;

addpath([code_directory, filesep, 'toolbox']);
addpath([code_directory, filesep, 'lib']);
addpath([code_directory, filesep, 'spm12']);

cutfreq=0.1164;  % cutoff frequency

% Load Neuromorphometric template
mask_mni = spm_read_vols(spm_vol([code_directory, filesep, 'labels_Neuromorphometrics', filesep, 'wlabels_Neuromorphometrics_unique.nii']));

mask_mni_sub = mask_mni;
mask_list_eff = unique(mask_mni);
mask_list_eff = mask_list_eff(2:end);
mask_mni_len = length(mask_list_eff);

% read header info
cd([wkdir, filesep, subname]);
para_file_ID = fopen('parameter_RS.txt', 'r');
while ~feof(para_file_ID)
    
    tline = fgetl(para_file_ID);
    if regexp(tline, 'SmoothFWHM')  %second line indicates the SmoothFWHmm        
        colon_loc =regexp(tline, ':');
        SmoothFWHMmm = str2num(tline(colon_loc+1:end));
    end
    
    if regexp(tline, 'TRs')  %second line indicates the SmoothFWHmm
        colon_loc =regexp(tline, ':');
        tr = str2num(tline(colon_loc+1:end));
    end
    
    if regexp(tline, 'boldFileName_RS')  %second line indicates the SmoothFWHmm
        colon_loc =regexp(tline, '"');
        [bold_dir, bold_file, bold_file_ext] = fileparts(tline(colon_loc(3)+1:colon_loc(4)-1));
    end
end

disp(tr)
fs=1/tr;
fcutoff = cutfreq/(fs/2);  %Hz
[b,a]= butter(2,fcutoff);

%%%%%%rp file
rp_name = dir([wkdir, filesep, subname, filesep, bold_dir, filesep, 'rp*.txt']);
filerp = fopen([rp_name.folder, filesep, rp_name.name],'r');
rp_data = fscanf(filerp,'%f %f %f %f %f %f',[6,inf])';
fclose(filerp);

for ii = 1:size(rp_data, 2)
    rp_data(:, ii) = filtfilt(b, a, detrend(rp_data(:, ii), 1));
end

% Smooth the bold image. 
% get realigned and resliced bold scan (prefixed with 'r')
cd([wkdir, filesep, subname, filesep, bold_dir]);
P = cell(1,1);
P{1}   = spm_select('List', [wkdir, filesep, subname, filesep, bold_dir], 'img', ['^drwar*']);

% get the scan's data
V       = spm_vol(P);
V       = cat(1,V{:});

disp('smoothing');

for ii = 1:length(V)
    [pth,nam,ext] = fileparts(V(ii).fname);
    fnameIn       = fullfile(pth, [nam ext]);
    fname         = fullfile(pth, ['s' int2str(SmoothFWHMmm) nam ext]);
    spm_smooth(fnameIn, fname, SmoothFWHMmm);
end

% clear V for further use
clear V;

% Load smoothed image
smoothFiles = spm_vol(spm_select('List', [wkdir, filesep, subname, filesep, bold_dir], 'img', ['^s' int2str(SmoothFWHMmm) 'drwar*']));

dynnum = length(smoothFiles);
boldImg =  zeros(mni_dimension(1), mni_dimension(2), mni_dimension(3), dynnum);

for ii = 1:length(smoothFiles)
    boldImg(:, :, :, ii) = spm_read_vols(spm_vol(smoothFiles(ii).fname));
end

boldImg_2D = reshape(boldImg, mni_dimension(1)*mni_dimension(2)*mni_dimension(3), dynnum);

%% Calculate the CVR map
cerebellumImg = spm_read_vols(spm_vol([wkdir, filesep, subname, filesep, 'mask', filesep, 'brainMask_RS_cerebellum.nii']));
cerebellum_mask_loc = find(cerebellumImg == 1);

% whole brain mask and seed mask
maskImg = spm_read_vols(spm_vol([wkdir, filesep, subname, filesep, 'mask', filesep, 'brainMask_RS.nii']));
mask_loc = find(maskImg == 1);

mask_mni_sub = maskImg.*mask_mni_sub; %outside the brain = 0

% filter on the within-brain voxel, then convert it back to 4D images
boldImg_2D_eff = boldImg_2D(mask_loc, :);
boldImg_2D_eff_filt = zeros(size(boldImg_2D_eff));
for ii = 1:size(boldImg_2D_eff, 1)
    boldImg_2D_eff_filt(ii, :) = filtfilt(b,a,boldImg_2D_eff(ii, :));  % Filtering
end

boldImg_3D_filt = zeros(mni_dimension(1)*mni_dimension(2)*mni_dimension(3), 1);
boldImg_4D_filt = zeros(mni_dimension(1)*mni_dimension(2)*mni_dimension(3), dynnum);
for ii = 1:size(boldImg_2D_eff_filt, 2)
    boldImg_3D_filt(mask_loc) = boldImg_2D_eff_filt(:, ii);
    boldImg_4D_filt(:, ii) = boldImg_3D_filt;
end

% Using the averaged cerebellum signal as the surrogate for EtCO2.
bestEtCO2 = mean(boldImg_4D_filt(cerebellum_mask_loc, :), 1)';
clear boldImg_4D_filt;

% CVR_beta0 = Mean bold image
CVR_beta0 = zeros(size(boldImg, 1)*size(boldImg, 2)*size(boldImg, 3), 1);
CVR_beta0(mask_loc) = mean(boldImg_2D_eff_filt, 2);

% Calculate CVR_beta1
CVR_beta1 = zeros(size(boldImg, 1)*size(boldImg, 2)*size(boldImg, 3), 1);

regression_x = cat(2, bestEtCO2, rp_data, rp_data.^2, ones(length(rp_data), 1));
regrssion_x_psu = (regression_x'*regression_x)\regression_x';
regression_para = regrssion_x_psu*boldImg_2D_eff_filt';
CVR_beta1(mask_loc) = regression_para(1, :);

% Regress out the global signal
boldImg_2D_eff_filt = boldImg_2D_eff_filt-(regression_x(:, 1:end-1)*regression_para(1:end-1, :))';
boldImg_2D_eff_filt = zscore(boldImg_2D_eff_filt, 0, 2);

boldImg_3D_filt_res = zeros(mni_dimension(1)*mni_dimension(2)*mni_dimension(3), 1);
boldImg_4D_filt_res = zeros(mni_dimension(1)*mni_dimension(2)*mni_dimension(3), dynnum);
for ii = 1:size(boldImg_2D_eff_filt, 2)
    boldImg_3D_filt_res(mask_loc) = boldImg_2D_eff_filt(:, ii);
    boldImg_4D_filt_res(:, ii) = boldImg_3D_filt_res;
end

% Average BOLD signal in each ROI mask
boldImg_roi_1D_res =  zeros(mask_mni_len, dynnum);
for ii = 1:mask_mni_len
    boldImg_roi_2D = boldImg_4D_filt_res(mask_mni_sub==mask_list_eff(ii), :);
    boldImg_roi_1D_res(ii, :) = mean(boldImg_roi_2D, 1);
    clear boldImg_roi_2D;
end

% Calculate correlation maps
corrMap = zeros(mni_dimension(1)*mni_dimension(2)*mni_dimension(3), mask_mni_len);
for uu = 1:mask_mni_len
    bestEtCO2_roi_res = boldImg_roi_1D_res(uu, :)';
    bestEtCO2_roi_res = zscore(bestEtCO2_roi_res, 0, 1);
    
    regression_x_roi = cat(2, bestEtCO2_roi_res, ones(length(rp_data), 1));
    CVR_roi = zeros(size(boldImg, 1)*size(boldImg, 2)*size(boldImg, 3), 1);
    
    regrssion_x_roi_psu = (regression_x_roi'*regression_x_roi)\regression_x_roi';
    beta_roi_para = regrssion_x_roi_psu*boldImg_2D_eff_filt';
    
    CVR_roi(mask_loc) = beta_roi_para(1, :)';
    corrMap(:, uu) = CVR_roi;
    clear CVR_roi
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Combine CVR beta maps and corrMaps together
CVR_4D = reshape(cat(2, corrMap, CVR_beta0, CVR_beta1), mni_dimension(1), mni_dimension(2), mni_dimension(3), size(corrMap, 2)+2);

batFile = dir([wkdir, filesep, subname, filesep, 'CVR_voxelshift_etco2_cerebellum', filesep, '*_s8_co2delay.img']);
batImg = spm_read_vols(spm_vol([batFile.folder, filesep, batFile.name]));
CVR_4D = cat(4, CVR_4D, batImg);

for uu = 1:size(CVR_4D, 4)
    CVR_4D(:,:,:,uu) = rewrite3D(CVR_4D(:,:,:,uu), maskImg);
end

mkdir([wkdir, filesep, subname, filesep, 'DLRS_input']);
prefix_mri = [wkdir, filesep, subname, filesep, 'DLRS_input', filesep, 'DLRS_input_layer_'];
img3Dto2D(CVR_4D, prefix_mri, 1, 0, 0);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function img_f = rewrite3D(img, mask)

img(isnan(img))=0;

mask_loc = find(mask == 1);
mask_loc0 = find(mask == 0);
img(mask_loc0) = nan;
img_mean = mean(img(mask_loc));
img_std = std(img(mask_loc), 0, 1);
img_re = (img - img_mean)/img_std;
img_f = img_re;
img_f(mask_loc0) = 0;
end

