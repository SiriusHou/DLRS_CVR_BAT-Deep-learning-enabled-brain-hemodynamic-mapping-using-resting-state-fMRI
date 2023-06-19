function relativeCvrMap = zachCVR_generateRelativeCvrMap(cvrDir,cvrFileName,brainMask)
% takes a cvr map directory, file name, and a brain mask (3D matrix). Then
% finds the global average cvr value of the voxels inside the brain mask.
% Finally, creates a new cvr map that contains the relative cvr values
% (relativeCvr = originalCvr/globalAverageCvr)

cvrMap = spm_read_vols(spm_vol([cvrDir filesep cvrFileName '.img']));

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