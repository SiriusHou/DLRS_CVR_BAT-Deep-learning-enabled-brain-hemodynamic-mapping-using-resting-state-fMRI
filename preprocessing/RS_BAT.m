function RS_BAT(wkdir, subname, code_directory)

cd([wkdir, filesep, subname]);
mni_resolution = [2, 2, 2];
mni_type = 16;

addpath([code_directory, filesep, 'toolbox']);
addpath([code_directory, filesep, 'lib']);
addpath([code_directory, filesep, 'spm12']);


cutfreq=0.1164;  % cutoff frequency

%%%%%%read header info
para_file_ID = fopen('parameter_RS.txt', 'r');
while ~feof(para_file_ID)
    
    tline = fgetl(para_file_ID);
    if regexp(tline, 'boldFileName_RS')  %second line indicates the SmoothFWHmm
        colon_loc =regexp(tline, '"');
        [bold_file_dir, bold_file_name, bold_file_ext] = fileparts(tline(colon_loc(3)+1:colon_loc(4)-1));
    end
    
    if regexp(tline, 'TRs')  %second line indicates the SmoothFWHmm
        colon_loc =regexp(tline, ':');
        tr = str2num(tline(colon_loc+1:end));
    end
end

cd([wkdir, filesep, subname, filesep, bold_file_dir]);

disp(tr)
fs=1/tr;
fcutoff = cutfreq/(fs/2);  %Hz
[b,a]= butter(2,fcutoff);

%%%%%%rp file
rp_name = dir([wkdir, filesep, subname, filesep, bold_file_dir, filesep, 'rp*.txt']);
filerp = fopen([rp_name.folder, filesep, rp_name.name],'r');
rp_data = fscanf(filerp,'%f %f %f %f %f %f',[6,inf])';
fclose(filerp);


%bold image
bold_image_name = dir([wkdir, filesep, subname, filesep, bold_file_dir, filesep, 'war*.nii']);
boldImg = spm_read_vols(spm_vol([bold_image_name.folder, filesep, bold_image_name.name]));
dynnum = size(boldImg, 4);

%% detrend
mask = spm_read_vols(spm_vol([wkdir, filesep, subname, filesep, 'mask', filesep, 'brainMask_RS.nii']));
brainVox = find(mask == 1);
detrendImg_4D = zeros(size(boldImg));
mask_size = size(mask);

for vox = 1:length(brainVox)
    [row,col,sl] = ind2sub(mask_size, brainVox(vox));
    TS1 = squeeze(boldImg(row,col,sl,:));
    
    sig = detrend(TS1);
    meansig = mean(TS1-sig);
    TS2 =  sig + meansig;     % detrended signal
    detrendImg_4D(row,col,sl,:) = TS2;
end

% write filtered and detrended images
ftempname = [wkdir, filesep, subname, filesep, bold_file_dir, filesep, 'dr' bold_image_name.name];
write_hdrimg(detrendImg_4D, ftempname, mni_resolution, mni_type);

[ftemp_dir, ftemp_name, ftemp_ext] = fileparts(ftempname);
for ii = 1:size(detrendImg_4D, 4)
    if ii < 10
        write_hdrimg(detrendImg_4D(:, :, :, ii), [ftemp_dir, filesep, ftemp_name, '-00', num2str(ii), '-001.img'], mni_resolution, mni_type); 
    elseif ii < 100
        write_hdrimg(detrendImg_4D(:, :, :, ii), [ftemp_dir, filesep, ftemp_name, '-0', num2str(ii), '-001.img'], mni_resolution, mni_type); 
    else
        write_hdrimg(detrendImg_4D(:, :, :, ii), [ftemp_dir, filesep, ftemp_name, '-', num2str(ii), '-001.img'], mni_resolution, mni_type); 
    end
end

boldImg_2D = reshape(detrendImg_4D, size(boldImg, 1)*size(boldImg, 2)*size(boldImg, 3), size(boldImg, 4));

%% cerebellum mask
grayImg = spm_read_vols(spm_vol([wkdir, filesep, subname, filesep, 'mask', filesep, 'brainMask_RS_cerebellum.nii']));
graymask_loc = find(grayImg == 1);

%whole brain mask and seed mask
maskImg = spm_read_vols(spm_vol([wkdir, filesep, subname, filesep, 'mask', filesep, 'brainMask_RS.nii']));
boldImg_4D_eff = boldImg_2D(graymask_loc, :);

%filter
boldImg_4D_eff_filt = zeros(size(boldImg_4D_eff));
for ii = 1:size(boldImg_4D_eff, 1)
    boldImg_4D_eff_filt(ii, :) = filtfilt(b,a,boldImg_4D_eff(ii, :));  % Filtering
end

bestEtCO2 = nanmean(boldImg_4D_eff_filt, 1)';
nave = round(length(bestEtCO2)/4);
[Y,I]=sort(bestEtCO2,'descend');
EtCO2_min =mean(Y(end-nave:end));  % lowest 1/4 as baseline
EtCO2_mean=mean(Y);
EtCO2_max =mean(Y(1:nave));        % highest 1/4 for output

avgBOLD_fname = 'avgBOLD_cerebellum.txt';
EtCO2_t = tr*linspace(1, length(bestEtCO2), length(bestEtCO2))';
fileID_avgBOLD = fopen(avgBOLD_fname,'w');
for zz = 1:length(EtCO2_t)
    fprintf(fileID_avgBOLD,'%f,%f\n', EtCO2_t(zz), bestEtCO2(zz));
end
fclose(fileID_avgBOLD);

% Calcualte RS BAT map
cvrdir_v = [wkdir, filesep, subname, filesep, 'CVR_voxelshift_etco2_cerebellum'];
mkdir(cvrdir_v);

fixedDelayRange(1) = -9;
fixedDelayRange(2) = 9;
SmoothFWHMmm=8;
voxelwiseResult = CVR_mapping_voxelwise_GLM_RS([wkdir, filesep, subname, filesep, bold_file_dir], ['dr' bold_image_name.name], tr, SmoothFWHMmm, fixedDelayRange, maskImg, ...
    avgBOLD_fname, EtCO2_mean, EtCO2_min, cvrdir_v, rp_data);

close all;
end
