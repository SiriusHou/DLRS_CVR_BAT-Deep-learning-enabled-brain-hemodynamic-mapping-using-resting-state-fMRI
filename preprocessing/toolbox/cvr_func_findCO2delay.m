function [optDelay, optEtCO2, minres, coefs] =...
    cvr_func_findCO2delay(boldcurve, tr, etco2timecourse,...
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
dynnum  = size(boldcurve,1);

% for ii = 1:11
%     etco2curve  = etco2plot_interp08(etco2timecourse,delays(ii),tr);
%     etco2curve  = etco2curve(1:length(boldcurve));
%     temp        = corrcoef(boldcurve',etco2curve);
%     cc(ii)      = temp(2,1);
% end
bbs = [];
for ii = 1:samples
    etco2curve = cvr_func_interpTimecourse(etco2timecourse,delays(ii),dynnum,tr);   
    [bb, Yp, Xp, res] = cvr_func_glm(etco2curve,boldcurve(:,2));
    resid(ii)  = res;
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
    [~, resid_loc] = min(resid);
    delay_min = delays(resid_loc);
    
    % fitting method ----------------------------------------------------------
    warning('off','all')
    pp       = polyfit(delays',resid,4);
    warning('on','all')
%     [x1,f1]  = fminbnd(@(x)(pp(1)*x^4+pp(2)*x^3+pp(3)*x^2+pp(4)*x+pp(5)),delay_range(1),delay_range(2));
    [x1,f1]  = fminbnd(@(x)(pp(1)*x^4+pp(2)*x^3+pp(3)*x^2+pp(4)*x+pp(5)),max(delay_min-2, delay_range(1)), min(delay_min+2, delay_range(2)));

    optDelay = x1;
    minres   = f1;
    Vq       = polyval(pp,delays2);
end

% if findAbsoluteMin
%     optDelay = optDelay + delay_range(1);
% end

% % interpolating method ----------------------------------------------------
% Vq      = interp1(delays',resid,delays2','spline');
% [minres, ind] = min(Vq);
% optDelay = delays2(ind);
% % figure, contour(Vq);
% % figure, imshow(Vq);

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


