function output = CVR_mapping_voxelwise_GLM_RS(varargin)
% CVR_mapping_voxelwise_GLM(path_temp,name_cvr,TR,SmoothFWHMmm, delayrange, brainmask,etco2name,outpath_temp);

path_temp = varargin{1};
name_cvr = varargin{2};
TR = varargin{3};
SmoothFWHMmm = varargin{4};
delaypara = varargin{5};
brainmask = varargin{6};
etco2FilePath = varargin{7};
EtCO2_mean = varargin{8};
EtCO2_min = varargin{9};
if nargin>9
    outpath_temp = varargin{10};
end
rp_data = varargin{11};

% matlabpool('open',4);
% p = parpool(4);

minDelay=delaypara(1);
maxDelay=delaypara(2);


P = spm_select('FPList',path_temp,['^', name_cvr]); %Altered YUS Dec 26, 2007   // might need to change
V = spm_vol(P);
nVol = size(V,1);
matsize = V(1).dim;
voxsize = abs(V(1).mat(eye(4,3)==1))';

etco2timecourse = dlmread(etco2FilePath, ',', 0, 0);

mask=reshape(brainmask,[size(brainmask,1)*size(brainmask,2)*size(brainmask,3),1]);
maskvox=find(mask==1);
voxbeta_0001=zeros(size(maskvox,1),1);
voxbeta_0002=zeros(size(maskvox,1),1);
voxbeta_0003=zeros(size(maskvox,1),1);
fitCC=zeros(size(maskvox,1),1);
img=zeros(length(maskvox),nVol);

cutfreq=0.1164;  % cutoff frequency
fs=1/TR;
fcutoff = cutfreq/(fs/2);  %Hz
[b,a]= butter(2,fcutoff);
    
for ii = 1:size(rp_data, 2)   
    rp_data(:, ii) = filtfilt(b, a, detrend(rp_data(:, ii), 1));
end

for i=1:nVol
    img_temp=spm_read_vols(V(i));
    img1=reshape(img_temp,[size(mask,1)*size(mask,2)*size(mask,3),1]);
    img2=img1(maskvox);
    img(:,i)=img2;
end;

w = waitbar(0,'Running Voxelwise Shift...');

tenth = floor(length(maskvox)./10);
sections = {1:tenth,...
            1+tenth:2*tenth,...
            1+2*tenth:3*tenth,...
            1+3*tenth:4*tenth,...
            1+4*tenth:5*tenth,...
            1+5*tenth:6*tenth,...
            1+6*tenth:7*tenth,...
            1+7*tenth:8*tenth,...
            1+8*tenth:9*tenth,...
            1+9*tenth:length(maskvox),...
            };

p = parpool(12);
for wctr = 1:length(sections)
    waitbar(wctr./length(sections),w,'Running Voxelwise Shift...')
    parfor k=sections{wctr} %101190:101190%
        sig = squeeze(filtfilt(b, a, img(k,:)))'; %low pass filter
        sig = smooth(sig, 10, 'sgolay', 5);
        sig = [TR*(0:length(sig)-1)',sig];
        % Initial        
        optDelay = cvr_func_findCO2delay_RS(sig, TR, etco2timecourse, rp_data, [minDelay maxDelay],0,0,1);
        voxdelay(k)=optDelay;
        voxoptEtCO2 = cvr_func_interpTimecourse(etco2timecourse,optDelay,nVol,TR);
        [coefs, Yp, ~, ~] = cvr_func_glm(voxoptEtCO2,sig(:,2));
        voxbeta_0001(k) = coefs(1);
        voxbeta_0002(k) = coefs(2);
        voxbeta_0003(k) = coefs(3);
        fitCC(k) = corr(Yp,sig(:,2));

    end
end

delete(p);
disp('Finished voxelwise shift.')
close(w)

hcdelayall=zeros(length(mask),1); 
voxbeta1all=zeros(length(mask),1);
voxbeta2all=zeros(length(mask),1);
voxbeta3all=zeros(length(mask),1);
fitccall=zeros(length(mask),1);

hcdelayall(maskvox)=voxdelay;
voxbeta1all(maskvox)=voxbeta_0001; 
voxbeta2all(maskvox)=voxbeta_0002; 
voxbeta3all(maskvox)=voxbeta_0003; 
fitccall(maskvox)=fitCC; 

hcdelaymap=reshape(hcdelayall,size(brainmask));
beta1map=reshape(voxbeta1all,size(brainmask));
beta2map=reshape(voxbeta2all,size(brainmask));
beta3map=reshape(voxbeta3all,size(brainmask));
CCmap=reshape(fitccall,size(brainmask));

bold_batPath = [outpath_temp filesep name_cvr '_co2delay.img'];
bold_beta1Path = [outpath_temp filesep name_cvr '_s' int2str(SmoothFWHMmm) '_beta_0001.img'];
bold_beta2Path = [outpath_temp filesep name_cvr '_s' int2str(SmoothFWHMmm) '_beta_0002.img'];
bold_beta3Path = [outpath_temp filesep name_cvr '_s' int2str(SmoothFWHMmm) '_beta_0003.img'];
bold_ccPath = [outpath_temp filesep name_cvr '_s' int2str(SmoothFWHMmm) '_cc.img'];

write_ANALYZE(hcdelaymap,bold_batPath,matsize,voxsize,1,16);
write_ANALYZE(beta1map,bold_beta1Path,matsize,voxsize,1,16);
write_ANALYZE(beta2map,bold_beta2Path,matsize,voxsize,1,16);
write_ANALYZE(beta3map,bold_beta3Path,matsize,voxsize,1,16);
write_ANALYZE(CCmap,bold_ccPath,matsize,voxsize,1,16);

spm_smooth(bold_batPath , [outpath_temp filesep name_cvr '_s' int2str(SmoothFWHMmm) '_co2delay.img'], SmoothFWHMmm);
% thre = 4;
% nave=floor(nVol/thre);
% avedelay=mean(voxdelay);
% optEtCO2 = cvr_func_interpTimecourse(etco2timecourse,avedelay,nVol,TR);
% [Yco2,I]=sort(optEtCO2,'descend');
% EtCO2_min =mean(Yco2(end-nave:end));  % lowest 1/4 as baseline
% EtCO2_mean=mean(Yco2);

co2_CVR=(beta1map./(beta3map-beta1map*(EtCO2_mean-EtCO2_min)))*100.*brainmask;
write_ANALYZE(co2_CVR,[outpath_temp filesep name_cvr '_s' int2str(SmoothFWHMmm) '_CVR.img'],matsize,voxsize,1,16);

% matlabpool('close');
% delete(p);

output = {bold_batPath, bold_beta1Path, bold_beta2Path, bold_beta3Path,...
    bold_ccPath, EtCO2_mean, EtCO2_min};

end


