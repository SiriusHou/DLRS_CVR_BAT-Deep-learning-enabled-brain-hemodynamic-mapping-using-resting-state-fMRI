function [optDelayEtCO2, optEtCO2, optDelayEtO2, optEtO2, minres, coefs] = ...
    cvr_func_findDelay_hcho(boldcurve,tr, etco2timecourse, eto2timecourse, delay_range,idraw)

% Use this script to find the delay for which the extracted etco2 
% has the best correlation to the average BOLD signal extracted from
% thalamus of normalised time series

stepdly = 0.02; % unite: s
samples = 11;
dynnum  = length(boldcurve);
delays  = linspace(delay_range(1),delay_range(2),samples);
delays2 = delay_range(1):stepdly:delay_range(2);
cc      = zeros(samples,samples);

all_etco2curve = zeros(dynnum,samples);
all_eto2curve = zeros(dynnum,samples);
for ii = 1:samples
    etco2tmp  = cvr_func_interpTimecourse(etco2timecourse,delays(ii),dynnum,tr); % use interp script to do the interp at each selected i for selected file  
    all_etco2curve(:,ii)  = etco2tmp(1:dynnum); % select length of interp etco2 to equal the no. of BOLD dyn

    eto2tmp  = cvr_func_interpTimecourse(eto2timecourse,delays(ii),dynnum,tr);
    all_eto2curve(:,ii)  = eto2tmp(1:dynnum);
end

for ii = 1:samples
for jj = 1:samples
    Y = boldcurve;
    X = [all_etco2curve(:,ii) all_eto2curve(:,jj)];
    [coefs, Yp, Xcomp, res] = cvr_func_glm(X,Y);
    cc(ii,jj) = res;
end
end
% figure, imshow(all_eto2curve);
% figure, imshow(cc);
% figure, contour(cc);

% interpolate to find the max cc
[X,Y]   = meshgrid(delays);
[Xq,Yq] = meshgrid(delays2);
V       = cc;
Vq      = interp2(X,Y,V,Xq,Yq,'spline');
% figure, contour(Vq,10); hold on; plot(ind_etco2,ind_eto2,'r*');
% figure, imshow(Vq); 

% [maxcc, ccind] = max(Vq(:));
[minres, ccind] = min(Vq(:));
ind_eto2 = floor((ccind - 1)/size(Vq,1)) + 1;
ind_etco2  = mod(int32(ccind - 1),size(Vq,2));
% Vq(ind_eto2,ind_etco2) = 0;

optDelayEtCO2 = delays2(ind_etco2);
optDelayEtO2  = delays2(ind_eto2);
optEtCO2  = cvr_func_interpTimecourse(etco2timecourse,optDelayEtCO2,dynnum,tr);
optEtO2   = cvr_func_interpTimecourse(eto2timecourse, optDelayEtO2, dynnum,tr);
optEtCO2  = optEtCO2(1:dynnum);
optEtO2   = optEtO2(1:dynnum);

% plot figures
[coefs, Yp, Xcomp, ~] = cvr_func_glm([optEtCO2 optEtO2],boldcurve);

if idraw == 1
    figure,
    subplot(1,2,1), contour(Vq), hold on, plot(ind_eto2,ind_etco2,'r*'); hold off;
    subplot(1,2,2), plot(1:dynnum, boldcurve, 1:dynnum, Yp,'r-','linewidth',2);
    
    % check fitting
    Ycom = Xcomp(:,1:2) + 0.9*[Xcomp(:,4) Xcomp(:,4)];
    figure, plot(1:dynnum, Yp,'Linewidth',2); hold on, plot(Ycom); hold off;
end


