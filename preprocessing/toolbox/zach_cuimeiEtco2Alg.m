function output = zach_cuimeiEtco2Alg(timeData,co2Data)
% takes a column of timestamps and a column of co2 levels and finds the
% end-tidal co2 peaks. This is done by smoothing the co2 trace and then
% using the crosses between the smoothed trace and the actual trace to
% search for local maxima.

    % get CO2t and COt out and get rid of the0
    CO2t=timeData;
    CO2=co2Data;
%     for idx=2:size(CO2t,1)
%         if (CO2t(idx)<0.001)
%             CO2t=CO2t(1:idx-1);
%             CO2=CO2(1:idx-1);
%             break;
%         end
%     end
    SmoothCO2=smooth(CO2,100); %10 is the smooth span need to otimize
    
    %looking for cross
    imax=0;
    for idx=2:size(CO2,1)-1
        if ((SmoothCO2(idx)-CO2(idx))>0)&((SmoothCO2(idx-1)-CO2(idx-1))<0)
            imax=imax+1;
            areaidx(imax)=idx;
            if (imax==1)

                 [peakCO2(imax,2), maxidx]=max(CO2(1:idx));
                 tmpCOt=CO2t(1:idx);
                  peakCO2(imax,1)=tmpCOt(maxidx);
            else
                 [peakCO2(imax,2), maxidx] =max(CO2(areaidx(imax-1):idx));
                    tmpCOt=CO2t(areaidx(imax-1):idx);
                  peakCO2(imax,1)=tmpCOt(maxidx);
            end
        end
    end
    
    output = peakCO2;

end