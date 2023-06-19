function result = zach_CVR_calc_betamaps(outPath, beta1Path, beta3Path,...
    EtCO2_mean, EtCO2_min, brainmaskPath)

% load data for caluclation
beta1 = loadimage(beta1Path);
beta3 = loadimage(beta3Path);
brnmsk = loadimage(brainmaskPath);

% calculate CVR
cvrMap=((beta1./(beta3-beta1*(EtCO2_mean-EtCO2_min)))*100).*brnmsk;

% write to .hdr/.img files cvrMap
p_out1  = spm_vol(beta1Path);
p_out1.fname = [outPath];
cvrMapPath = p_out1.fname;
spm_write_vol(p_out1,cvrMap);

relMap = zachCVR_generateRelativeCvrMap(cvrMap,brnmsk);

% write to .hdr/.img files for relMap
p_out2  = spm_vol(beta1Path);
[pth,nme,ext] = fileparts(outPath);
p_out2.fname = [pth filesep 'rel_' nme ext];
relCvrMapPath = p_out2.fname;
spm_write_vol(p_out2,relMap);

% give results
result = {cvrMapPath,relCvrMapPath};


function relativeCvrMap = zachCVR_generateRelativeCvrMap(cvrMap, brainMask)
% finds the global average cvr value of the voxels inside the brain mask
% and creates a new cvr map containing the relative values
% (relativeCvr = originalCvr/globalAverageCvr)

voxelCtr = 0;
total = 0;
for i = 1:size(cvrMap,1)
    for ii = 1:size(cvrMap,2)
        for iii = 1:size(cvrMap,3)
            if(brainMask(i,ii,iii))
                total = total+cvrMap(i,ii,iii);
                voxelCtr = voxelCtr + 1;
            end
        end
    end
end

averageCvr = total/voxelCtr;

relativeCvrMap = cvrMap./averageCvr;

end

end