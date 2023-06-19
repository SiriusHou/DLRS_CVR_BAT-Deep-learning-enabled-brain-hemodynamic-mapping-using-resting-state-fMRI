function [matrixSize,matrixCenter,res] = zach_get_matSize_res_center(imgPath,imgFile)
%%
% takes the location of a .img or .hdr file and returns the scan's avergae 
% 3D matrix size and center, as well as the average 3D resolution in mm


% imgPath = 'C:\Users\Zach\OneDrive\College - Grad\Lu Lab\peiyingCvrCode\workspace';
% imgFile = 'hc_bold.img';
%%

    % gets scan information
    V = spm_vol([imgPath filesep imgFile]);

    % gets affline transformation matricies and dimensions of each scan
    transMats = {V.mat};
    dimVecs = {V.dim};

    % finds average affline transformation matrix and average dimension
    % (dimensions should be consistent throughout scan)
    totalTrans = zeros(4);
    totalDims = zeros(1,3);
    for i = 1:length(transMats)
        totalTrans = totalTrans + transMats{i};
        totalDims = totalDims + dimVecs{i};
    end
    avgTransMat = totalTrans./i;
    avgDimVec = totalDims./i;
    
    % gets matrixSize and matrixCenter
    matrixSize = avgDimVec;
    matrixCenter = matrixSize./2;
    
    % gets resolution in mm from the affline transformation matrix;
    % resolution is indicated by the diagonal (negative values and the
    % time component are omitted)
    res = abs([avgTransMat(1,1),avgTransMat(2,2),avgTransMat(3,3)]);

% end