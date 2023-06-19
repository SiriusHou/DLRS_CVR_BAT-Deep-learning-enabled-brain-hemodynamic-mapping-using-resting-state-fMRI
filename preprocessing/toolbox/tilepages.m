function plane = tilepages(vols, varargin)
% Show 4D image arrays across dynamics with multislices
% Show 3D image arrays across slices on a row*col mat
% yli, 20141014

[rr,cc,pp,vv] = size(vols);

if nargin > 1
    Nr = varargin{1};
    if isempty(varargin{2})
        Nc = ceil(pp/Nr);
    else
        Nc = varargin{2};
    end
else
    Nr = floor(sqrt(pp));
    Nc = ceil(pp/Nr);
end

if vv == 1 % 3D image
    tmp   = reshape(vols,rr,cc*pp);
    tmp   = [tmp zeros(rr,Nr*Nc*cc-cc*pp)];
    tmp   = reshape(tmp,rr,Nc*cc,Nr);
    tmp   = rotpages90(tmp,1);
    tmp   = reshape(tmp,Nc*cc,rr*Nr);
    plane = rot90(tmp,3);
    % figure, imshow(tmp);
else % 4D images
    tmp   = rotpages90(vols,1);
    tmp   = reshape(tmp,cc,rr*pp,vv);
    tmp   = rotpages90( tmp,3);
    plane = reshape(tmp,rr*pp,vv*cc);
    % figure, imshow(plane);
end

