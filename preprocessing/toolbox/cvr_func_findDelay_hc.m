function [optDelay, optEtCO2, minres, coefs] = cvr_func_findDelay_hc(boldtimecourse, etco2timecourse, delay_range,idraw)
% Use this script to find the delay for which the extracted etco2 
% has the best correlation to the average BOLD signal extracted from
% thalamus of normalised time series the difference is, this algorithm will
% not do constant filling if the interpolation exceeds the range

%% setup
samples = 16; % how many delays will be considered
stepdly = 0.02; % Unit: 0.5 s
delays  = linspace(delay_range(1),delay_range(2),samples);
resid   = zeros(samples,1); % stores residuals

% store values associated with (average?) BOLD timecourse
boldcurve = boldtimecourse(:,2);
timeline  = boldtimecourse(:,1);
dynnum  = length(boldcurve);
tr      = mean(timeline(2:end) - timeline(1:end-1));

% for each delay, find residual value
for ii = 1:samples
    % format etco2timecourse according to delay
    etco2curve = cvr_func_interpTimecourse(etco2timecourse,delays(ii),dynnum,tr); 
    % get residual value (res)
    [bb, Yp, Xp, res] = cvr_func_glm(etco2curve,boldcurve);
    % save res value
    resid(ii)  = res;
end

% fit to find the min res
delays2 = linspace(delay_range(1),delay_range(2),abs(delay_range(2)-delay_range(1))/stepdly);

% fitting method ----------------------------------------------------------
% 4th degree polynomial fitted to the residual curve (you'll probably find
    % it to resemble a parabola)
pp       = polyfit(delays',resid,4);
[x1,f1]  = fminbnd(@(x)(pp(1)*x^4+pp(2)*x^3+pp(3)*x^2+pp(4)*x+pp(5)),...
    delay_range(1),delay_range(2));
optDelay = x1;
minres   = f1;
Vq       = polyval(pp,delays2);

% get ouputs
optEtCO2 = cvr_func_interpTimecourse(etco2timecourse,optDelay,dynnum,tr);
[coefs, Yp, ~, ~] = cvr_func_glm(optEtCO2,boldcurve);

% draw data if specified in input by idraw
if idraw == 1
    origEtCO2 = cvr_func_interpTimecourse(etco2timecourse,0,dynnum,tr);
    X = origEtCO2;
    [rr, ~] = size(X);    
    X1 = [X, (1:rr)']; % add linear terms
    X2 = X1 - repmat(mean(X1,1),rr,1);
    X2 = [X2, ones(rr,1)];
    Yori = X2 * coefs;

    figure('pos',[300 300 1000 300]),
    subplot(1,3,    1), plot(delays',resid,'k+',delays2',Vq,'k-',optDelay,minres,'r*');
                    ylabel('Fitting residual'); xlabel('EtCO_2 time course shift');
                    title('Fitting residual minimum found');
    subplot(1,3,[2 3]), plot(boldcurve,'Linewidth',2); hold on, 
                    plot(Yori,'-','Linewidth',2,'color',[1 1 1]*0.8); 
                    plot(Yp,'r-','Linewidth',2); hold off;
                    title('Fitting residual minimum found');
                    legend('BOLD','EtCO2 (original)','EtCO2 (shifted)');
                    xlim([0 dynnum]); ylabel('BOLD signal'); xlabel('Dyn #');
                    title('BOLD time course fitted by EtCO_2 time course');
    tightfig;
end


