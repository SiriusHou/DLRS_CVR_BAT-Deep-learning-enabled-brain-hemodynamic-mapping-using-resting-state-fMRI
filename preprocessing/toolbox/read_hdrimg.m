function [imgvols,varargout] = read_hdrimg(imgfname)
% read and return hdr/img files

vs = spm_vol(imgfname);
imgvols = spm_read_vols(vs); % vs.mat dont have to be same
imgvols(isnan(imgvols)) = 0; % remove nan

matsize = vs(1).dim;
voxsize = abs(vs(1).mat([1,6,11]));
dt      = vs(1).dt(1);
ss      = vs(1).pinfo(1);

if nargout == 1
    
elseif nargout == 2
    varargout{1} = vs;
elseif nargout == 3
    varargout{1} = matsize;
    varargout{2} = voxsize;
elseif nargout == 5
    varargout{1} = matsize;
    varargout{2} = voxsize;
    varargout{3} = dt;
    varargout{4} = ss;
else
    disp('Number of outputs is not as designed...');
end

