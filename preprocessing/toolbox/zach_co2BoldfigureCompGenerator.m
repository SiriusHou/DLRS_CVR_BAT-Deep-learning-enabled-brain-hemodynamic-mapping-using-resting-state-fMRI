function zach_co2BoldfigureCompGenerator(co2Path,boldPath,outDir,TR)
    
%     co2Path = spm_select('FPList',[inPath filesep 'workspace'],'boldSynced.*$');
%     boldPath = spm_select('FPList',[inPath filesep 'workspace'],'WB_ns.*$');
    

    
    fid1 = fopen(co2Path);
    fid2 = fopen(boldPath);
    co2Data = textscan(fid1,'%f%f','delimiter',',');
    co2Data = [co2Data{1},co2Data{2}];
    boldData = textscan(fid2,'%f');
    boldData = boldData{1};
    boldData = [TR*(0:length(boldData)-1)',boldData];
    
    fclose(fid1)
    fclose(fid2)
    
%     co2Data(:,2) = co2Data(:,2)-mean(co2Data(:,2));
%     boldData(:,2) = boldData(:,2)-mean(boldData(:,2));
    
%     [~,nme] = fileparts(inPath);
    
    f = figure;
    [hAx,hLine1,hLine2] = plotyy(boldData(:,1),boldData(:,2),co2Data(:,1),co2Data(:,2));
    xlabel('Time (s)')
    ylabel(hAx(1),'BOLD')
    ylabel(hAx(2),'EtCO2 (mmHg)')
    set(hLine1,'linewidth',2)
    set(hLine2,'linewidth',2)
    set(hLine1,'color','b')
    set(hLine2,'color','r')
    set(hAx(1),'ycolor','b')
    set(hAx(2),'ycolor','r')

    grid on
    legend('BOLD','EtCO2')

%     saveas(f,[outPath filesep nme '_etco2Bold_compare.fig'])
    saveas(f,[outDir filesep 'etco2Bold_compare.fig'])
    
end
    