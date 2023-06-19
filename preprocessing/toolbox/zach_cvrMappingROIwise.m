function  zach_cvrMappingROIwise(mniBoldDir, mniMaskFile, etco2Dir, tr,...
    delayRange,t1Path,outDir)

    % need:
    % mni bold scans | need mniBoldDir to point to the MNI dir in dumpDir
    % mni brain mask | mnimskfile is the path to the 2mm resolution mni msk
    % etco2 | path filename and ext
    % tr
    % delay range
    % outDir
    % t1Path and name_T1

%     % NO DOWNSAMPLE VVVVV
%     
%     % get mni sections from t1 data and downsample 289 label mni image
%     P0 = spm_select('FPList',t1Path,'_289Labels.*_MNI.img$');
%     if isempty(P0) || size(P0,1) ~= 1
%         P0 = [t1Path filesep name_T1 '_mr_286Labels_MNI.img'];
%     end
%     P1 = [P0(1:end-4) '_2mm.img'];
%     ASL_downSampleMNI(P0,P1);
%     
%     % ^^^^^^
    
    
    % extract roi values
    bold_19roi = [];
    bold_54roi = [];
    bold_289roi = [];
    
    bold_mni_files = spm_select('FPList',mniBoldDir,'MNI.img$');
    
    rrbold = bold_mni_files;

    [sum1,sum2,sum3,name1,name2,name3,nvox1,nvox2,nvox3,...
        missingVox1,missingVox2,missingVox3] =...
        cvr_t1roi_boldtimecourse_mni(t1Path,rrbold,mniMaskFile);
    bold_19roi = [bold_19roi, sum1];
    bold_54roi = [bold_54roi, sum2];
    bold_289roi = [bold_289roi, sum3];
    
    ratioIn1 = nvox1./(nvox1+missingVox1);
    ratioIn2 = nvox2./(nvox2+missingVox2);
    ratioIn3 = nvox3./(nvox3+missingVox3);
    
    % export roi cvr results
    fresult = [outDir filesep 'cvr_roi_results.txt'];
    fid = fopen(fresult,'wt');
    coltitle = 'CVR results by ROIs (Number of ROI = %d)\n';
    colnames = 'Index\tROI_name\tNumber_of_voxels\tRatio_of_roi_included\tCVR(%%BOLD/mmHg)\tCO2_BAT(sec)\tBOLD_EtCO2_cc\tNotes\n';
    colformt = '%d\t%s\t%u\t%1.2f\t%1.4f\t%1.2f\t%1.4f\t%s\n';
    
    fprintf(fid,coltitle,19);
    fprintf(fid,colnames);
    for iroi = 1:19
        if iroi < size(bold_19roi,1) + 1
            roisig = bold_19roi(iroi,:)';
            [roi_cvr,roi_delay,roi_cc] = CVR_mapping_roi_GLM(roisig,etco2Dir,tr,delayRange);
            if roi_cc <= 0.6
                warnMsg = 'LOW CC';
            else
                warnMsg = '';
            end
            if ratioIn1(iroi)<0.5
                warnMsg = ['Insf Ratio ' warnMsg];
                fprintf(fid,colformt,iroi,name1{iroi},nvox1(iroi),ratioIn1(iroi),NaN,NaN,NaN,warnMsg);
            else
                fprintf(fid,colformt,iroi,name1{iroi},nvox1(iroi),ratioIn1(iroi),roi_cvr,roi_delay,roi_cc,warnMsg);
            end
        end
    end
    fprintf(fid,'\n\n');
    fprintf(fid,coltitle,54);
    fprintf(fid,colnames);
    for iroi = 1:54
        if iroi < size(bold_54roi,1) + 1
            roisig = bold_54roi(iroi,:)';
            [roi_cvr,roi_delay,roi_cc] = CVR_mapping_roi_GLM(roisig,etco2Dir,tr,delayRange);
            if roi_cc <= 0.6
                warnMsg = 'LOW CC';
            else
                warnMsg = '';
            end
            if ratioIn2(iroi)<0.5
                warnMsg = ['Insf Ratio ' warnMsg];
                fprintf(fid,colformt,iroi,name2{iroi},nvox2(iroi),ratioIn2(iroi),NaN,NaN,NaN,warnMsg);
            else
                fprintf(fid,colformt,iroi,name2{iroi},nvox2(iroi),ratioIn2(iroi),roi_cvr,roi_delay,roi_cc,warnMsg);
            end
        end
    end
    fprintf(fid,'\n\n');
    fprintf(fid,coltitle,289);
    fprintf(fid,colnames);
    for iroi = 1:289
        if iroi < size(bold_289roi,1) + 1
            roisig = bold_289roi(iroi,:)';
            [roi_cvr,roi_delay,roi_cc] = CVR_mapping_roi_GLM(roisig,etco2Dir,tr,delayRange);
            if roi_cc <= 0.6
                warnMsg = 'LOW CC';
            else
                warnMsg = '';
            end
            if ratioIn3(iroi)<0.5
                warnMsg = ['Insf Ratio ' warnMsg];
                fprintf(fid,colformt,iroi,name3{iroi},nvox3(iroi),ratioIn3(iroi),NaN,NaN,NaN,warnMsg);
            else
                fprintf(fid,colformt,iroi,name3{iroi},nvox3(iroi),ratioIn3(iroi),roi_cvr,roi_delay,roi_cc,warnMsg);
            end
        end
    end
    fprintf(fid,'\n\n');
    fclose(fid);
    
end