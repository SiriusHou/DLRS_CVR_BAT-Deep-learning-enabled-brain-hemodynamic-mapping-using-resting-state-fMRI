function mniDir = zach_mniNormalize3DFileDir(dumpDir, commandLocation, t1Path)

    % takes location (path) of a directory containing 3D .img and .hdr
    % files corresponding to a 4D brain scan
    % creates a new set of .img and .hdr files in MNI space, in a new
    % directory in dumpDir titled MNI
    
    % Requires the directory containing the two files needed
    % to get the MNI normalization (IMG_apply_AIR_tform1 and
    % IMG_change_res_info)


    outDir = [dumpDir filesep 'MNI'];
    mkdir(outDir)
    scans = spm_select('FPList',dumpDir,'ts\d*.img$');
    
    for i = 1:length(scans)
    
        scanDir = scans(i,:);
        
        [scanPath, scanName] = fileparts(scanDir);
    
        outPath = [outDir filesep scanName '_MNI.img'];
%         outPath2mm = [outDir filesep scanName '_MNI2mm.img'];
    
        % sepcify terminal commands
        cmd1 = [commandLocation filesep 'IMG_apply_AIR_tform1 ' ...
                scanDir ' ' outPath ' ' ...
                t1Path filesep 'matrix_air.txt 1 ' ...
                t1Path filesep 'mni.imgsize 1'];
        cmd2 = [commandLocation filesep 'IMG_change_res_info ' ...
                outPath ' ' outPath ' 1 1 1'];

        system(cmd1)
        system(cmd2)

%         ASL_downSampleMNI(outPath,outPath2mm)
    
    end
    
    mniDir = outDir;

end