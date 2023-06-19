function [optDelay, optEtCO2, minres, coefs] =...
    cvr_func_findCO2delay_RS(boldcurve, tr, etco2timecourse, rp_data, ...
    delay_range,idraw,findAbsoluteMin,delayIteration,varargin)

% Use this script to find the delay for which the extracted etco2 
% has the best correlation to the average BOLD signal extracted from
% thalamus of normalised time series the difference is, this algorithm will
% not do constant filling if the interpolation exceeds the range


if ~isempty(varargin)
    outDirFigures = varargin{1};
end

if length(varargin) > 1
    globalBoldFlag = varargin{2};
else
    globalBoldFlag = 0;
end

samples = length(delay_range(1):delayIteration:delay_range(2));
% samples = 31;
stepdly = 0.01; % Unit: 0.5 s
delays  = linspace(delay_range(1),delay_range(2),samples);
resid   = zeros(samples,1);
index   = zeros(samples,1);
% index1   = zeros(samples,1);

coeff_co2  = zeros(samples,1);
dynnum  = size(boldcurve,1);

bbs = [];
for ii = 1:samples
    etco2curve = cvr_func_interpTimecourse(etco2timecourse,delays(ii),dynnum,tr); 
    [bb, Yp, Xp, res, se] = cvr_func_glm_RS(etco2curve, boldcurve(:,2), rp_data);
   
    coeff_co2(ii) = bb(1)/se(1);
    
    resid(ii)  = res;
    index(ii) = coeff_co2(ii);

    if coeff_co2(ii)<0
        resid(ii) = nan;
    end
    bbs = [bbs,bb];
end
% plot(bbs)

if findAbsoluteMin
    [minres,optDelayIdx] = min(resid);
    % etco2's coefficient should not be negative when finding global BOLD
    % shift. If this is the case, select the next-lowest local minimum
    if bbs(1,optDelayIdx) < 0 && globalBoldFlag
        [pks,locs] = findpeaks(-resid);
        pks = -pks;
        [srtPks,I] = sort(pks);
        srtLocs = locs(I);
        check = 0;
        for i = 1:length(srtLocs)
            if bbs(1,srtLocs(i)) > 0
                check = 1;
                minres = srtPks(i);
                optDelayIdx = srtLocs(i);
                break
            end
        end
        if ~check
            error('No local minima found in regression. Check global BOLD and EtCO2 signal quality.')
        end
    end
    optDelay = delays(optDelayIdx);
else
    % fit to find the max cc
    delays2 = delay_range(1) : stepdly : delay_range(2);
    if length(find(index<=0)) == length(index)
        resid_loc = length(index);
    else
        [~, resid_loc] = max(index,[],'omitnan');
    end

%     [~, resid_loc] = min(resid);
    delay_min = delays(resid_loc);
    
    % fitting method ----------------------------------------------------------
    warning('off','all')
    resid_loc = find(~isnan(resid));  
    pp       = polyfit(delays(resid_loc)',resid(resid_loc),4);
    warning('on','all')
       
    [x1,f1]  = fminbnd(@(x)(pp(1)*x^4+pp(2)*x^3+pp(3)*x^2+pp(4)*x+pp(5)),...
    max(delay_min-delayIteration, delay_range(1)), min(delay_min+delayIteration, delay_range(2)));

    optDelay = x1;
    minres   = f1;
    Vq       = polyval(pp,delays2);
end

optEtCO2 = cvr_func_interpTimecourse(etco2timecourse,optDelay,dynnum,tr);
[coefs, Yp, ~, ~] = cvr_func_glm(optEtCO2,boldcurve(:,2));


if idraw == 1
    zach_showShift(boldcurve,etco2timecourse,[optDelay+etco2timecourse(:,1),etco2timecourse(:,2)])
    f = figure;
    if findAbsoluteMin
        subplot(1,2,1), plot(delays',resid,'ko',optDelay,minres,'r*');
    else
        subplot(1,2,1), plot(delays',resid,'ko',delays2',Vq,'.-',optDelay,minres,'r*');
    end
    subplot(1,2,2), plot(boldcurve(:,2),'Linewidth',2); hold on, plot(Yp,'r-','Linewidth',2); hold off;
    if findAbsoluteMin
        saveas(f,[outDirFigures filesep 'boldEtco2Residuals_shift1.fig'])
    else
        saveas(f,[outDirFigures filesep 'boldEtco2Residuals_shift2.fig'])
    end
end


