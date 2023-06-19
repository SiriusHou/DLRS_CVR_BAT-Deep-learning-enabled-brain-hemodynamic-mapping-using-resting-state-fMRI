function [outVol,varargout] = reorientVol(inVol,ddd,varargin)
% This function recorients image array
% from arbitrary order to +x+y+z order:
%
% dim  dir  cords       dim  dir  cords
% -------------------------------------
%  1    -     z          2    +     x
%  2    +     x   --->   3    +     y
%  3    -     y          1    +     z
%
% or from +x+y+z order to arbitrary order if inverse:
%
% dim  dir  cords       dim  dir  cords
% -------------------------------------
%  1    +     x          3    -     z
%  2    +     y   --->   1    +     x
%  3    +     z          2    -     y
%
% Input:
%    inVol : 3D matrix containing image intensity
%    ddd : 6-digit character showing the direction of 1st, 2nd, and
%    3rd dimension (e.g. '+y-z-x')
%    l_inverse : if -1, perform the inverse conversion
%
% Output:
%    outVol : three-dimensional matrix of the image where the 1st,
%              2nd and 3rd dimensions goes x, y, and z directions
%    mat_xyz : 1x3 integer matrix, the sizes of x, y, and z dimensions.
%    order_xyz : 1x3 integer matrix, the dimension of the original image
%                matrix corresponding to x, y, and z respectively.
%    sign_xyz : 1x3 integer matrix, the sign of the original image array
%               direction w.r.t. x, y, and z.
%
% e.g.
% In the Philip Achieva system, the x, y, z directions goes (right-hand coordinates):
%  1st dim : -x -> +x : R -> L
%  2nd dim : -y -> +y : A -> P
%  3rd dim : -z -> +z : F -> H
%
% But in Analyze format (left-hand coordinates),
%  1st dim : -x -> +x : R -> L
%  2nd dim : -y -> +y : P -> A (!)
%  3rd dim : -z -> +z : F -> H
%
% My coronal image from Achieva has the following dimension directions:
% 1st dim : R -> L
% 2nd dim : H -> F
% 3rd dim : A -> P
% 
% If you want a radiological view in Analyze format,
% [OutArray,mat_xyz,order_xyz,sign_xyz] =reorientVol(InArray,'+x-z-y')
% The output outArray gives the image matrix where each dimension
% goes from R to L (+x), from P to A (+y), and from F to H (+z).
% mat_xyz = [mat1 mat3 mat2], order_xyz = [1 3 2], sign_xyz = [1 -1 -1]
%
% yli20160715

if length(ddd) == 6
    ddd = reshape(ddd,2,3)';
else
    disp('Wrong input...');
end

if nargin == 1
    l_inverse = varargin{1};
else
    l_inverse = 1;
end

% tic
if l_inverse == 1
    tmp1 = [['1';'2';'3'] ddd];
elseif l_inverse == -1
    tmp1 = ['1+x','2+y','3+z'];
    tmp1 = strrep(tmp1,['+' ddd(1,2)],[ddd(1,1) 'a']);
    tmp1 = strrep(tmp1,['+' ddd(2,2)],[ddd(2,1) 'b']);
    tmp1 = strrep(tmp1,['+' ddd(3,2)],[ddd(3,1) 'c']);
    tmp1 = reshape(tmp1,3,3)';
end

tmp2 = sortrows(tmp1,3);
xyz_order = str2num(tmp2(:,1))';
xyz_sign  = (tmp2(:,2)=='+')';

outVol = inVol;
outVol = permute(outVol,xyz_order);
for ss = 1:3
    if xyz_sign(ss)==0
        outVol = flip(outVol,ss);
    end
end
% toc

varargout{1} = size(outVol);
varargout{2} = xyz_order;
varargout{3} = xyz_sign;

