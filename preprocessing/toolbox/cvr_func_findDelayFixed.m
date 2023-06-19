function [optDelay, optEtCO2, optEtO2, minres, coefs] =...
    cvr_func_findDelayFixed(boldcurve, tr, etco2timecourse, eto2timecourse,...
    delay_range,idraw,findAbsoluteMin,delayIteration,varargin)

% Use this script to find the delay for which the extracted etco2 
% has the best correlation to the average BOLD signal extracted from
% thalamus of normalised time series the difference is, this algorithm will
% not do constant filling if the interpolation exceeds the range


if ~isempty(varargin)
    outDirFigures = varargin{1};
end

samples = length(delay_range(1):delayIteration:delay_range(2));
delays  = linspace(delay_range(1),delay_range(2),samples);

% samples_o2 = length(delay_range_o2(1):delayIteration:delay_range_o2(2));
% delays_o2  = linspace(delay_range_o2(1),delay_range_o2(2),samples_o2);

resid = zeros(samples,1);
stepdly = 0.01;
dynnum  = length(boldcurve);

for ii = 1:samples
    etco2curve = cvr_func_interpTimecourse(etco2timecourse,delays(ii),dynnum,tr);
    eto2curve = cvr_func_interpTimecourse(eto2timecourse,delays(ii),dynnum,tr);
    [bb, Yp, Xp, res] = zach_cvr_func_glm_xyz(etco2curve,eto2curve,...
        boldcurve);
    resid(ii) = res;
end


if findAbsoluteMin
    
    [minres,optDelayIdx] = min(resid);
    optDelay = delays(optDelayIdx);
    
%     [minRow,minCol] = zach_getMinCoords(resid);
%     optDelayCo2 = delays_co2(minRow);
%     optDelayO2 = delays_o2(minCol);
%     minres = resid(minRow,minCol);
    
else
    % fit to find the max cc
    [~, resid_loc] = min(resid);
    delay_min = delays(resid_loc);
    delays2 = delay_range(1) : stepdly : delay_range(2);

    % fitting method ----------------------------------------------------------
    warning('off','all')
    pp       = polyfit(delays',resid,4);
    warning('on','all')
    [x1,f1]  = fminbnd(@(x)(pp(1)*x^4+pp(2)*x^3+pp(3)*x^2+...
        pp(4)*x+pp(5)),max(delay_min-2, delay_range(1)), min(delay_min+2, delay_range(2)));
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
optEtO2 = cvr_func_interpTimecourse(eto2timecourse,optDelay,dynnum,tr);
[coefs, Yp, ~, ~] = zach_cvr_func_glm_xyz(optEtCO2,optEtO2,boldcurve);

if idraw == 1
    zach_showShift2([(tr*(0:length(boldcurve)-1))',boldcurve],...
        etco2timecourse,[optDelay+etco2timecourse(:,1),...
        etco2timecourse(:,2)],eto2timecourse,...
        [optDelay+eto2timecourse(:,1),eto2timecourse(:,2)])
%     legend('BOLD','CO_2','','O_2','')
    f = figure;
    if findAbsoluteMin
        subplot(1,2,1), plot(delays',resid,'ko',optDelay,minres,'r*');
    else
        subplot(1,2,1), plot(delays',resid,'ko',delays2',Vq,'.-',optDelay,minres,'r*');
    end
    subplot(1,2,2), plot(boldcurve,'Linewidth',2); hold on, plot(Yp,'r-','Linewidth',2); hold off;
%     if findAbsoluteMin
%         saveas(f,[outDirFigures filesep 'boldEtco2Residuals_shift1.fig'])
%     else
%         saveas(f,[outDirFigures filesep 'boldEtco2Residuals_shift2.fig'])
%     end
end

% if idraw == 1
%     zach_showShift2([(tr*(0:length(boldcurve)-1))',boldcurve],...
%         etco2timecourse,[optDelay+etco2timecourse(:,1),...
%         etco2timecourse(:,2)],eto2timecourse,...
%         [optDelay+eto2timecourse(:,1),eto2timecourse(:,2)])
%     legend('BOLD','CO_2','','O_2','')
%     f = figure;
%     if findAbsoluteMin
%         [X,Y] = meshgrid(delays,delays);
%         subplot(1,2,1); surf(X,Y,resid);
%         xlabel('CO_2 Delay');ylabel('O_2 Delay');zlabel('Residuals');
%         hold on
%         scatter3(optDelay,optDelayCo2,minres,'r*','linewidth',5);
% %         plot(delays_co2',resid,'ko',optDelayCo2,minres,'r*');
%     else
%         subplot(1,2,1), plot(delays_co2',resid,'ko',delays2',Vq,'.-',optDelayCo2,minres,'r*');
%     end
%     subplot(1,2,2), plot(boldcurve,'Linewidth',2); hold on, plot(Yp,'r-','Linewidth',2); hold off;
% %     subplot(2,2,2)%, plot(boldcurve,'Linewidth',2); hold on, plot(Yp,'r-','Linewidth',2); hold off;
% %     subplot(2,2,4)%,
%     if findAbsoluteMin
%         saveas(f,[outDirFigures filesep 'boldEtco2Residuals_shift1.fig'])
%     else
%         saveas(f,[outDirFigures filesep 'boldEtco2Residuals_shift2.fig'])
%     end
% end


