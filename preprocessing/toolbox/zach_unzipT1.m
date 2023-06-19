function [t1Path,name_T1] = zach_unzipT1(workDir,multiatlasFileName)

    unzipPath = [workDir filesep multiatlasFileName(1:end-4)];

    if exist(unzipPath,'dir')
        rmdir(unzipPath,'s')
    end
    mkdir(unzipPath)
    
    unzip([unzipPath '.zip'], unzipPath)
    
    % find lowest point
    loc = unzipPath;
    name = multiatlasFileName(1:end-4);
    while 1
        
        x = dir(loc);
        folders = 0;
        for i = 3:length(x)
            if x(i).isdir
                folders = folders+1;
                loc = [loc filesep x(i).name];
                name = x(i).name;
            end
        end
        if folders > 1
            error('T1 Segmentation zip may not have subdirectories with multiple folders')
        end
        if folders == 0
            break
        end
        
        
    end
    
    t1Path = loc;
    name_T1 = name;




end