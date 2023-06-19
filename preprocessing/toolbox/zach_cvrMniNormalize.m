function zach_cvrMniNormalize(mainPath,cvrDir,path_out,cvrName,cvrOut_a,cvrOut_r,t1Path)
    % normalizes CVR maps (absolute and relative) to MNI space using cmd
    % comands

    disp('ASLMRICloud: Normalize CVR maps to MNI space...');
    
    % set up temporary paths and filenames for the MNI maps
    tmpCvrMni_a = [cvrDir filesep cvrName '_aCVR_MNI.img'];
    tmpCvrMni_r = [cvrDir filesep cvrName '_rCVR_MNI.img'];
    
    % set up output paths and filenames for the MNI maps
    outCvrMni_a = [path_out filesep cvrName '_aCVR_MNI.img'];
    outCvrMni_r = [path_out filesep cvrName '_rCVR_MNI.img'];
    
    % sepcify terminal commands for the absolute cvr map
    cmd1_a = [mainPath filesep 'toolbox' filesep 'IMG_apply_AIR_tform1 ' ...
            cvrOut_a ' ' tmpCvrMni_a ' ' ...
            t1Path filesep 'matrix_air.txt 1 ' ...
            t1Path filesep 'mni.imgsize 1'];
    cmd2_a = [mainPath filesep 'toolbox' filesep 'IMG_change_res_info ' ...
            tmpCvrMni_a ' ' tmpCvrMni_a ' 1 1 1'];
    
    % specify terminal commands for the relative cvr map
    cmd1_r = [mainPath filesep 'toolbox' filesep 'IMG_apply_AIR_tform1 ' ...
        cvrOut_r ' ' tmpCvrMni_r ' ' ...
        t1Path filesep 'matrix_air.txt 1 ' ...
        t1Path filesep 'mni.imgsize 1'];
    cmd2_r = [mainPath filesep 'toolbox' filesep 'IMG_change_res_info ' ...
            tmpCvrMni_r ' ' tmpCvrMni_r ' 1 1 1'];
    
    % run terminal commands for the absolute cvr map, then downsample
    system(cmd1_a);
    system(cmd2_a);
    ASL_downSampleMNI(tmpCvrMni_a,outCvrMni_a);
    
    % run terminal commands for the relative cvr map, then downsample
    system(cmd1_r);
    system(cmd2_r);
    ASL_downSampleMNI(tmpCvrMni_r,outCvrMni_r);
    
    
end
