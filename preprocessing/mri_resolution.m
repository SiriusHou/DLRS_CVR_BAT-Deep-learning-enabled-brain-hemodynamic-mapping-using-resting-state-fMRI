function resolution = mri_resolution(M)
    M = M(1:3, 1:3);
    resolution = sqrt(sum(M.^2, 1));
end