function paths = zach_split4DScan(pathTo4DScan,dumpDir)
    % takes a path to a 4D scan and produces n .img and .hdr files, where n is
    % the number of timepoints in the 4D scan. Places files in dumpDir
    % This functuion is dependent on spm12.
    % returns cell of strings containing all paths of the 3D scans

    [~,fileName] = fileparts(pathTo4DScan);
    
    disp(['Splitting ' fileName ' into 3D scans'])
    
    V = spm_vol(pathTo4DScan);
    data = spm_read_vols(V);
    paths = cell(length(V),1);
    
    disp('...')

    places = length(num2str(length(V)));
    
    for i = 1:length(V)
        num = num2str(i);
        extraZeros = places-length(num);
        for j = 1:extraZeros
            num = ['0' num];
        end
        V(i).fname = [dumpDir filesep fileName '_ts' num '.img'];
        paths{i} = V(i).fname;
        V(i).n = [1 1];
        spm_write_vol(V(i),data(:,:,:,i));
    end
    
end