function [coefs, Yp, Xcomp, res, se] = cvr_func_glm_RS(X,Y,Z)

warning off

[rr, cc] = size(X);
[r1, c1] = size(Y);

Z_diff = cat(1, zeros(1,size(Z, 2)), diff(Z));
X_lin = linspace(1, length(X), length(X))';

Z = Z - mean(Z, 1);
X1 = [X, Z, Z.^2]; % add linear and quadratic terms

X2 = X1 - repmat(mean(X1,1),rr,1);
X3 = [X2, ones(rr,1)]; % add const terms

bb = X3 \ Y; % estimate coefs

Yp = X3 * bb; % estimated Y

if c1 == 1 % if only one voxel
    Xcomp = X3 .* repmat(bb',rr,1); % estimate contribution of each variables
else
    Xcomp = 0;
end

coefs = bb;

dY = Y - Yp; % calculate fitting residual
res = sum( dY(:).^2 );

cov = res/(rr-size(X3, 2))*inv(X3'*X3);
se = sqrt(diag(cov));

return



