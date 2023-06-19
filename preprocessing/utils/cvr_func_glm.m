function [coefs, Yp, Xcomp, res] = cvr_func_glm(X,Y)
% solving b = X \ Y, where X is added by linear & const terms
% output: coefs, coefs
%         Yp, estimated Y
%         Xp, estimated components from each col of X
%         res, fitting residual
%
% yli199@jhmi.edu
% 20170712

[rr, cc] = size(X);
[r1, c1] = size(Y);

X1 = [X, (1:rr)']; % add linear terms
X2 = X1 - repmat(mean(X1,1),rr,1);
X2 = [X2, ones(rr,1)]; % add const terms

bb = X2 \ Y; % estimate coefs

Yp = X2 * bb; % estimated Y

if c1 == 1 % if only one voxel
    Xcomp = X2 .* repmat(bb',rr,1); % estimate contribution of each variables
else
    Xcomp = 0;
end

coefs = bb;

dY = Y - Yp; % calculate fitting residual
res = sum( dY(:).^2 );

