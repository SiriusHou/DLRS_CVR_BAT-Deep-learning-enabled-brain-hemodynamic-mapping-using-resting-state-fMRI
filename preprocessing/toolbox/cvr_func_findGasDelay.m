function [optDelayCo2, optEtCO2, optDelayO2, optEtO2, minres, coefs] =...
    cvr_func_findGasDelay(boldcurve, tr, etco2timecourse, eto2timecourse,...
    delay_range_co2,delay_range_o2,idraw,findAbsoluteMin,delayIteration,varargin)

% Use this script to find the delay for which the extracted etco2 
% has the best correlation to the average BOLD signal extracted from
% thalamus of normalised time series the difference is, this algorithm will
% not do constant filling if the interpolation exceeds the range


if ~isempty(varargin)
    outDirFigures = varargin{1};
end

samples_co2 = length(delay_range_co2(1):delayIteration:delay_range_co2(2));
delays_co2  = linspace(delay_range_co2(1),delay_range_co2(2),samples_co2);

samples_o2 = length(delay_range_o2(1):delayIteration:delay_range_o2(2));
delays_o2  = linspace(delay_range_o2(1),delay_range_o2(2),samples_o2);

resid = zeros(samples_co2,samples_o2);
stepdly = 0.01;
dynnum  = length(boldcurve);

resid = resid;
w = waitbar(0);
for ii = 1:samples_co2
    waitbar(ii/samples_co2)
    for jj = 1:samples_o2
        etco2curve = cvr_func_interpTimecourse(etco2timecourse,delays_co2(ii),dynnum,tr);
        eto2curve = cvr_func_interpTimecourse(eto2timecourse,delays_o2(jj),dynnum,tr);
        [bb, Yp, Xp, res] = zach_cvr_func_glm_xyz(etco2curve,eto2curve,...
            boldcurve);
        resid(ii,jj)  = res;
    end
end
close(w)

[X,Y] = meshgrid(delays_o2,delays_co2);

if findAbsoluteMin
    [minRow,minCol] = zach_getMinCoords(resid);
    optDelayCo2 = delays_co2(minRow);
    optDelayO2 = delays_o2(minCol);
    minres = resid(minRow,minCol);
    
else
    Z = reshape(resid,[],1);
    XY = [reshape(X,[],1),reshape(Y,[],1)];
    F = fit(XY,Z,'poly22');
    x0 = [delays_o2(length(round(delays_o2/2))),...
        delays_co2(length(round(delays_co2/2)))];
%   x0 = [0,0];
    
    myFunc = (@(x,y) F.p00 + F.p10.*x + F.p01.*y + F.p20.*x.^2 +...
        F.p11.*x.*y + F.p02.*y.^2);
%     myFunc = (@(x,y) F.p00 + F.p10.*x + F.p01.*y + F.p20*x.^2 +...
%         F.p11.*x.*y + F.p02*y.^2 + F.p30.*x.^3 + F.p21.*x.^2.*y + F.p12.*x.*y.^2 +...
%         F.p03.*y.^3 + F.p40.*x.^4 + F.p31.*x.^3.*y + F.p22.*x.^2.*y.^2 +...
%         F.p13.*x.*y.^3 + F.p04.*y.^4);

    [x1,f1] = fminsearch(@(x0) myFunc(x0(1),x0(2)),x0);
    
    optDelayO2 = x1(1);
    optDelayCo2 = x1(2);
    minres   = f1;
    
%     % fit to find the max cc
%     delays2 = delay_range_co2(1) : stepdly : delay_range_co2(2);
% 
%     % fitting method ----------------------------------------------------------
%     warning('off','all')
%     pp       = polyfit(delays_co2',resid,4);
%     warning('on','all')
%     [x1,f1]  = fminbnd(@(x)(pp(1)*x^4+pp(2)*x^3+pp(3)*x^2+pp(4)*x+pp(5)),delay_range_co2(1),delay_range_co2(2));
%     optDelayCo2 = x1;
%     minres   = f1;
%     Vq       = polyval(pp,delays2);
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

optEtCO2 = cvr_func_interpTimecourse(etco2timecourse,optDelayCo2,dynnum,tr);
optEtO2 = cvr_func_interpTimecourse(eto2timecourse,optDelayO2,dynnum,tr);
[coefs, Yp, ~, ~] = zach_cvr_func_glm_xyz(optEtCO2,optEtO2,boldcurve);


if idraw == 1
    f = zach_showShift2([(tr*(0:length(boldcurve)-1))',boldcurve],...
        etco2timecourse,[optDelayCo2+etco2timecourse(:,1),...
        etco2timecourse(:,2)],eto2timecourse,...
        [optDelayO2+eto2timecourse(:,1),eto2timecourse(:,2)]);
    legend('BOLD','CO_2','','O_2','')
    f = figure;
    if findAbsoluteMin
        subplot(1,2,1); surf(X,Y,resid);
%         xlabel('O_2 Delay');ylabel('CO_2 Delay');zlabel('Residuals');
%         hold on
%         scatter3(optDelayO2,optDelayCo2,minres,'r*','linewidth',5);
%         plot(delays_co2',resid,'ko',optDelayCo2,minres,'r*');
    else
        dim = size(X);
        xx = reshape(X,[],1);
        yy = reshape(Y,[],1);
        funOut = myFunc(xx,yy);
        FUNOUT = reshape(funOut,dim(1),dim(2));
        
        subplot(1,2,1); surf(X,Y,FUNOUT);
%         subplot(1,2,1); scatter3(xxx,yyy,funOut);

%         scatter3()
%         subplot(1,2,1), plot(delays_co2',resid,'ko',delays2',Vq,'.-',optDelayCo2,minres,'r*');
    end
    xlabel('O_2 Delay');ylabel('CO_2 Delay');zlabel('Residuals');
    hold on
    scatter3(optDelayO2,optDelayCo2,minres,'r*','linewidth',5);
    subplot(1,2,2), plot(boldcurve,'Linewidth',2); hold on, plot(Yp,'r-','Linewidth',2); hold off;
%     subplot(2,2,2)%, plot(boldcurve,'Linewidth',2); hold on, plot(Yp,'r-','Linewidth',2); hold off;
%     subplot(2,2,4)%,
    if findAbsoluteMin
        saveas(f,[outDirFigures filesep 'boldEtco2Residuals_shift1.fig'])
    else
        saveas(f,[outDirFigures filesep 'boldEtco2Residuals_shift2.fig'])
    end
end


