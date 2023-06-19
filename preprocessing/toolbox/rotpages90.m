function outvols = rotpages90(vols,num)
% Rotate each page num*90 degree in a vols
% yli, 20141014
if nargin == 1
    num = 1;
end

[m1, m2, m3, m4] = size(vols);
outvols       = zeros(m2, m1, m3, m4);

for ivol = 1:m4
for islc = 1:m3
    slctmp = squeeze(vols(:,:,islc,ivol));
	outvols(:,:,islc,ivol) = rot90(slctmp,num);
end
end

