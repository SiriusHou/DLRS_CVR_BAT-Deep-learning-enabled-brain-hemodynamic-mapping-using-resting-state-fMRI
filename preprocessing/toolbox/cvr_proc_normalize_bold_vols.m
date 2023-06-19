% normalize bold vols to mpr to mni
% bold*281 -> bold_mpr*281 -> bold_mni*281

% clear variables;

for isub = 1:length(subid)
    subdir  = [studypath filesep 'hlu_vcid' subid{isub}]
    subname = spm_select('List',subdir,['^hlu_vcid' subid{isub} '_cvr-d0.*.img$']);
    subname = subname(1:end-4);
    
    % manage folder and directory
    subdirtmp = [subdir filesep 'tmp'];
%     subdirslt = [subdir filesep 'cvr_result'];
    mkdir(subdirtmp);
%     mkdir(subdirslt);
    
    % get mprage name from T1 result
    path_mpr    = [subdir filesep 'mpr'];
    name_tmp    = spm_select('List',path_mpr,'.*[^mni].imgsize$');
    name_mpr    = name_tmp(1:end-8);

    % split bold series to single vol files
    boldname = [subdir filesep 'r' subname '.img'];
    p_bold = spm_vol(boldname);
    v_bold = spm_read_vols(p_bold);
    matsize = size(v_bold);
    other   = cell(matsize(4),1);
    
    source = [subdir filesep 'mean' subname '.img'];
    p_src = spm_vol(source); % make sure source and other have same .mat
    
    for ivol = 1:matsize(4)
        p_out  = p_bold(1);
        p_out.fname = [subdirtmp filesep 'r' subname '-ivol' sprintf('%03d',ivol) '.img'];
        p_out.mat   = p_src.mat; % make sure source and other have same .mat
        
        spm_write_vol(p_out,v_bold(:,:,:,ivol));
        
        other{ivol,1} = p_out.fname;
    end
    
    % coregister all vols to mpr
%     mpr_brain = ASL_mprageSkullstrip(path_mpr, name_mpr);
%     target = mpr_brain;
    target = spm_select('FPList',path_mpr,'.*_brain.img$');
    source = [subdir filesep 'mean' subname '.img'];
    asl_coreg12(target,source,other);
    
    % normalize rbold to mni
    path_code = ['/home/yli/Documents/MATLAB'];
    rmaps = spm_select('FPList',subdirtmp,'^rrhlu_vcid.*-ivol.*[^mm].img$');
    for ii = 1:matsize(4)
        rmap = rmaps(ii,:);
        rmap_mni = [rmap(1:end-4) '_mni.img'];
        rmap_mni2= [rmap(1:end-4) '_mni_2mm.img'];
        
        cmd1 = [path_code filesep 'IMG_apply_AIR_tform1 ' ...
            rmap ' ' rmap_mni ' ' ...
            path_mpr filesep 'matrix_air.txt 1 ' path_mpr filesep 'mni.imgsize 1'];
        cmd2 = [path_code filesep 'IMG_change_res_info ' ...
            rmap_mni ' ' rmap_mni ' 1 1 1'];
        
        system(cmd1);
        system(cmd2);
        ASL_downSampleMNI(rmap_mni,rmap_mni2);
    end

    % remove bulky intermediate files
    delete([subdirtmp filesep '*_mni.*']);
    delete([subdirtmp filesep 'rhlu*ivol*']);

end



