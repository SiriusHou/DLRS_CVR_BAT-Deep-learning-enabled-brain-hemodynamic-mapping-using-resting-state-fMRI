function outPath = zach_mniNormalize(scanDir, commandLocation, t1Path)

    % converts .img file to mni space with (1x1x1)mm resolution

    [scanPath, scanName] = fileparts(scanDir);
    outDir = scanPath;

    outPath = [outDir filesep scanName '_MNI.img'];
%     outPath2mm = [outDir filesep scanName '_MNI2mm.img'];

    % sepcify terminal commands
    cmd1 = [commandLocation filesep 'IMG_apply_AIR_tform1 ' ...
            scanDir ' ' outPath ' ' ...
            t1Path filesep 'matrix_air.txt 1 ' ...
            t1Path filesep 'mni.imgsize 1'];
    cmd2 = [commandLocation filesep 'IMG_change_res_info ' ...
            outPath ' ' outPath ' 1 1 1'];

    system(cmd1)
    system(cmd2)

    % no downsampling!
%     ASL_downSampleMNI(outPath,outPath2mm)

end