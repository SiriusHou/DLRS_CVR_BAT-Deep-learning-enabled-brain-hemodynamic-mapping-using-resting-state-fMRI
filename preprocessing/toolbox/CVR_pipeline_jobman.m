function success = CVR_pipeline_jobman(path_code, jsonInput)

% path_code = 'C:\Users\Zach\OneDrive\College_Grad\Lu_Lab\peiyingYangStitchedPipeline_v2\commands';
% jsonInput = 'C:\Users\Zach\OneDrive\College_Grad\Lu_Lab\peiyingYangStitchedPipeline_v2\cvr_parameters.json';

    cvr_paras = loadjson(jsonInput);
    
    path_data = cvr_paras.Dir.mainPath;
    
    outputFileName = 'Result';
    
    % create a .log file to record output info/error msg
    datetimestr     = sprintf('%4d%02d%02d_%02d%02d%02d',int16(clock));
    diaryFile       = [path_data,filesep,'CVRMRICloud_' datetimestr '.log'];
    if ~isempty(dir([path_data,filesep,'*.log']))
        delete([path_data,filesep,'*.log']);
    end
    diary(diaryFile); % start the diary
    
    path_output = cvr_paras.Dir.outPath;

    success = false;
    

    try % if error happens, provide a .log for user to download
    %% single-subject processing only (for now)

        disp('CVRMRICloud: Single dataset...');
        
        tic
        output_content = CVR_pipeline(path_code, jsonInput);
        toc
        
        % zip the folder for download
        delete(   [path_output,filesep,outputFileName,'.zip']);
        zipname = [path_output,filesep,outputFileName,'.zip'];
        zip(zipname,{'*.txt','*.hdr','*.img','*.fig'},output_content);
        disp('CVRMRICloud zipped for download...');
        
        success = true;
        
        
        
    catch ME % if error happens, provide a .log for user to download
%         rethrow(ME)
        msgErr = getReport(ME);
        fileID = fopen(diaryFile,'a');
        fprintf(fileID,'\n');
        fprintf(fileID,'CVRMRICloud: Error message...\n');
        fprintf(fileID,'%s\n',msgErr);
        fclose(fileID);
        if ~isempty(dir([path_output,filesep,'*.log'])) % delete .log in path_output if any before copy
            delete([path_output,filesep,'*.log']);
        end
        copyfile(diaryFile, path_output);

        % zip .log for download
        zipname = [path_output,filesep,outputFileName,'.zip'];
        zip(zipname,diaryFile);
    end

diary off; % end the diary

