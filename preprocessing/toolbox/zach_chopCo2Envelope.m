function envelopeCo2FilePath = zach_chopCo2Envelope(co2FilePath,...
    avgBoldFilePath, envelope_interp_rate, TR)
    
    % takes co2 envelope, average bold signal, co2 recorder sample rate,
    % and repetition time to chop the co2 envelope to match with the 


    fid = fopen(co2FilePath);
    co2Data = textscan(fid,'%f%f','delimiter',',');
    co2Data = [co2Data{1},co2Data{2}];
    co2envlp = co2Data(:,2);
    
    fclose(fid);
    
    fid = fopen(avgBoldFilePath);
    bold = textscan(fid,'%f');
    bold = bold{1};
    fclose(fid);
    

    co2Time = (0:length(co2envlp)-1)';
    co2Time = co2Time./envelope_interp_rate;

    boldTime = (0:length(bold)-1)';
    boldTime = boldTime.*TR;
    boldTimeInterp = (0:1/envelope_interp_rate:max(boldTime))';

    boldInterp = interp1(boldTime,bold,boldTimeInterp);
    boldInterp(isnan(boldInterp))=[];

    figure;
    plot(boldTimeInterp,boldInterp-mean(bold),'r')
    hold on
    plot(co2Time,co2envlp-mean(co2envlp),'b')
    xlabel('Time (s)')

    
    
    lag = zach_subSignalXCorr(co2envlp-mean(co2envlp),boldInterp-mean(boldInterp));
    
%     [r,lags] = xcorr(co2envlp-mean(co2envlp),boldInterp-mean(boldInterp));
%     [~,idx] = max(r);
%     lag = lags(idx);
%     
    
    
    
    % if lag is neg then the user may not have sufficient data in the
    % co2 recording
    if lag < 0
%         envelopeCo2FilePath   = co2FilePath;
%         return
        error(['Error, see ' mfilename])
    end
    
    choppedCo2 = co2Data(lag:lag+length(boldTimeInterp),:);
    choppedCo2(:,1) = choppedCo2(:,1)-choppedCo2(1,1);
    figure;
    plot(choppedCo2(:,1),choppedCo2(:,2))
    title('chopped Co_2 Envelope')
    
    dlmwrite([co2FilePath(1:end-4) 'chopped.txt'], choppedCo2);

    envelopeCo2FilePath   = [co2FilePath(1:end-4) 'chopped.txt'];

    
    function lag = zach_subSignalXCorr(signal,subSignal)
        
        % finds the cross correlation, but only where a window the size of
        % the subsignal fits inside the signal. Then determines the lag
        % associated with the maximum cross correlation coefficient
        
        if length(subSignal) > length(signal)
            error('subSignal must be smaller than signal')
        end
        lag = 0;
        maxCoef = -inf;
        for i = 1:length(signal)-length(subSignal)+1
            
            windowedSignal = signal(i:i+length(subSignal)-1);
            
            x1 = windowedSignal;
            x2 = subSignal;
                        
            coef = (1/length(x1)) * sum(((x1-mean(x1))./std(x1)).*((x2-mean(x2))./std(x2)));
            
            if coef > maxCoef
                maxCoef = coef;
                lag = i-1;
            end
        end

    end
    
end