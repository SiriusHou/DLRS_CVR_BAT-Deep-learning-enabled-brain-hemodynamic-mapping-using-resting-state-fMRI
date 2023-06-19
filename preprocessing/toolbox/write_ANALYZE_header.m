
% This function writes only a header in ANALYZE format
%
% Input:
%    filename : The name of the file including the directory
%    mat_123    : 1 x 3 integer, the size of the matrix [mat1,mat2,mat3]
%    resolution : 1 x 3 integer, resolution of the image
%                 The default value is [1 1 1]
%
%    data_type: 4 for 'int16', 2 for 'uint8', 16 for 'float', 132 'uint16'
%               default is 4
%    mat_origin: 1 x 3 integer, the origin of the matrix
%               default is round(mat_123/2)
%
% Output:
%    The header file is stored
%    
% written by Jinsoo Uh (jinsoo.uh@utsouthwestern.edu)
% 2007-04-16

function f = write_ANALYZE_header(varargin)

filename = varargin{1};
mat_123_tmp = varargin{2};
mat_123 = reshape(mat_123_tmp,1,3);
resolution = varargin{3};
if nargin > 3
   data_type = varargin{4};
else
   data_type = 4;
end;
if nargin > 4
   mat_origin = varargin{5};
else
   mat_origin = round(mat_123/2);
end;

[fname_path fname_body fname_ext] = fileparts(filename);

    P = fullfile(fname_path, [fname_body '.hdr']);

    DIM = mat_123;
    VOX = resolution;
    SCALE = 1;
    TYPE = data_type;
    OFFSET = 0;
    ORIGIN = mat_origin;
    DESCRIP = '';
    write_analyze75_header(P,DIM,VOX,SCALE,TYPE,OFFSET,ORIGIN,DESCRIP);
        
    
