function [coefs, Yp, Xcomp, res, se] = zach_cvr_func_glm_xyz_XH(X,Y,Z)
% solving b = X \ Y, where X is added by linear & const terms
% output: coefs, coefs
%         Yp, estimated Y
%         Xp, estimated components from each col of X
%         res, fitting residual
%         se, add the standard error associated with each coefficient by Xirui Hou at 20201203
% yli199@jhmi.edu
% 20170712

[rr, cc] = size(X);
[r1, c1] = size(Z);

X1 = [X, Y, (1:rr)']; % add linear terms
X2 = X1 - repmat(mean(X1,1),rr,1);
X3 = [X2, ones(rr,1)]; % add const terms

bb = X3 \ Z; % estimate coefs

Yp = X3 * bb; % estimated Y

if c1 == 1 % if only one voxel
    Xcomp = X3 .* repmat(bb',rr,1); % estimate contribution of each variables
else
    Xcomp = 0;
end

coefs = bb;

dY = Z - Yp; % calculate fitting residual
res = sum( dY(:).^2 );

cov = res/(rr-size(X3, 2))*inv(X3'*X3);
se = sqrt(diag(cov));

return

