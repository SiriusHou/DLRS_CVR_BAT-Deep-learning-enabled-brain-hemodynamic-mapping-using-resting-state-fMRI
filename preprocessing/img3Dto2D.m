function img3Dto2D(CVR_g_4D, prefix_mri, d, ud, lr)

mni_resolution = [2, 2, 2];
mni_type = 16;

if d == 1 %z direction
    
    z_ini = 1;
    z_end = 91;
    
    for kk = z_ini:z_end
        
        bold_img = squeeze(CVR_g_4D(:,:, kk, :));
        
        if ud == 0 && lr == 0
            
            bold_img_f = padarray(bold_img,[2 1],0,'pre');
            bold_img_f = padarray(bold_img_f,[3 2],0,'post');
            
            if kk<10
                delete([prefix_mri, '00', num2str(kk), '.nii']);
                write_hdrimg(bold_img_f, [prefix_mri, '00', num2str(kk), '.nii'], mni_resolution, mni_type);
            else
                delete([prefix_mri, '0', num2str(kk), '.nii']);
                write_hdrimg(bold_img_f, [prefix_mri, '0', num2str(kk), '.nii'], mni_resolution, mni_type);
            end
        end
        
        if ud == 1 && lr == 0
            bold_img = flipud(bold_img);
            bold_img_f = padarray(bold_img,[2 1],0,'pre');
            bold_img_f = padarray(bold_img_f,[3 2],0,'post');
            if kk<10
                delete([prefix_mri, '10', num2str(kk), '.nii']);
                write_hdrimg(bold_img_f, [prefix_mri, '10', num2str(kk), '.nii'], mni_resolution, mni_type);
            else
                delete([prefix_mri, '1', num2str(kk), '.nii']);
                write_hdrimg(bold_img_f, [prefix_mri, '1', num2str(kk), '.nii'], mni_resolution, mni_type);
            end
        end
        
        if ud == 0 && lr == 1
            bold_img = fliplr(bold_img);
            bold_img_f = padarray(bold_img,[2 1],0,'pre');
            bold_img_f = padarray(bold_img_f,[3 2],0,'post');
            if kk<10
                delete([prefix_mri, '20', num2str(kk), '.nii']);
                write_hdrimg(bold_img_f, [prefix_mri, '20', num2str(kk), '.nii'], mni_resolution, mni_type);
            else
                delete([prefix_mri, '2', num2str(kk), '.nii']);
                write_hdrimg(bold_img_f, [prefix_mri, '2', num2str(kk), '.nii'], mni_resolution, mni_type);
            end
        end
        
    end
end

if d == 2 %y direction (coronal write)
    
    y_ini = 15;
    y_end = 95;
    
    for kk = y_ini:y_end
        
        bold_img = squeeze(CVR_g_4D(:, kk, :, :));
        
        if ud == 0 && lr == 0
            
            bold_img_f = padarray(bold_img,[2 9+1],0,'pre');
            bold_img_f = padarray(bold_img_f,[3 9+2],0,'post');
            
            if kk<10
                delete([prefix_mri, '00', num2str(kk), 'y.nii']);
                write_hdrimg(bold_img_f, [prefix_mri, '00', num2str(kk), 'y.nii'], mni_resolution, mni_type);
            else
                delete([prefix_mri, '0', num2str(kk), 'y.nii']);
                write_hdrimg(bold_img_f, [prefix_mri, '0', num2str(kk), 'y.nii'], mni_resolution, mni_type);
            end
        end
        
        if ud == 1 && lr == 0
            bold_img = flipud(bold_img);
            bold_img_f = padarray(bold_img,[2 9+1],0,'pre');
            bold_img_f = padarray(bold_img_f,[3 9+2],0,'post');
            if kk<10
                delete([prefix_mri, '10', num2str(kk), 'y.nii']);
                write_hdrimg(bold_img_f, [prefix_mri, '10', num2str(kk), 'y.nii'], mni_resolution, mni_type);
            else
                delete([prefix_mri, '1', num2str(kk), 'y.nii']);
                write_hdrimg(bold_img_f, [prefix_mri, '1', num2str(kk), 'y.nii'], mni_resolution, mni_type);
            end
        end
        
        if ud == 0 && lr == 1
            bold_img = fliplr(bold_img);
            bold_img_f = padarray(bold_img,[2 9+1],0,'pre');
            bold_img_f = padarray(bold_img_f,[3 9+2],0,'post');
            if kk<10
                delete([prefix_mri, '20', num2str(kk), 'y.nii']);
                write_hdrimg(bold_img_f, [prefix_mri, '20', num2str(kk), 'y.nii'], mni_resolution, mni_type);
            else
                delete([prefix_mri, '2', num2str(kk), 'y.nii']);
                write_hdrimg(bold_img_f, [prefix_mri, '2', num2str(kk), 'y.nii'], mni_resolution, mni_type);
            end
        end
        
    end
end


if d == 3 %x direction (sagittal write)
    
    x_ini = 15;
    x_end = 74;
    
    for kk = x_ini:x_end
        
        bold_img = permute(squeeze(CVR_g_4D(kk, :, :, :)), [2 1 3]);
        
        if ud == 0 && lr == 0
            
            bold_img_f = padarray(bold_img,[2 1],0,'pre');
            bold_img_f = padarray(bold_img_f,[3 2],0,'post');
            
            if kk<10
                delete([prefix_mri, '00', num2str(kk), 'z.nii']);
                write_hdrimg(bold_img_f, [prefix_mri, '00', num2str(kk), 'z.nii'], mni_resolution, mni_type);
            else
                delete([prefix_mri, '0', num2str(kk), 'z.nii']);
                write_hdrimg(bold_img_f, [prefix_mri, '0', num2str(kk), 'z.nii'], mni_resolution, mni_type);
            end
        end
        
        if ud == 1 && lr == 0
            bold_img = flipud(bold_img);
            bold_img_f = padarray(bold_img,[2 1],0,'pre');
            bold_img_f = padarray(bold_img_f,[3 2],0,'post');
            if kk<10
                delete([prefix_mri, '10', num2str(kk), 'z.nii']);
                write_hdrimg(bold_img_f, [prefix_mri, '10', num2str(kk), 'z.nii'], mni_resolution, mni_type);
            else
                delete([prefix_mri, '1', num2str(kk), 'z.nii']);
                write_hdrimg(bold_img_f, [prefix_mri, '1', num2str(kk), 'z.nii'], mni_resolution, mni_type);
            end
        end
        
        if ud == 0 && lr == 1
            bold_img = fliplr(bold_img);
            bold_img_f = padarray(bold_img,[2 1],0,'pre');
            bold_img_f = padarray(bold_img_f,[3 2],0,'post');
            if kk<10
                delete([prefix_mri, '20', num2str(kk), 'z.nii']);
                write_hdrimg(bold_img_f, [prefix_mri, '20', num2str(kk), 'z.nii'], mni_resolution, mni_type);
            else
                delete([prefix_mri, '2', num2str(kk), 'z.nii']);
                write_hdrimg(bold_img_f, [prefix_mri, '2', num2str(kk), 'z.nii'], mni_resolution, mni_type);
            end
        end
        
    end
end

end
