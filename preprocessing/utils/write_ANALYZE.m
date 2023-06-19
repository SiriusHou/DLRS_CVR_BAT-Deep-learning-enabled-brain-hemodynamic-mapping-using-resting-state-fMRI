
% This function writes a 3D image array into ANALYZE file with 'int16'
% format. 
%
% Input:
%    ImageArray : 3D matrix, the image array
%    filename : The name of the file.
%    mat_123    : 1 x 3 integer, the size of the matrix [mat1,mat2,mat3]
%    resolution : 1 x 3 float, resolution of the image
%                 The default value is [1 1 1]
%    scale_factor : The image array is scaled by this number before saved
%                   This is often useful when the image array has float 
%                   numbers and they lose significant digits if converted
%                   to integers.
%                   The default value is 1.
%    data_type : 4 ~ int16, 2 ~ uint8, 16 ~ float, default is 4
%    no_header:  if this is 1, the header file is not written
%                default is 0
%    mat_origin: 1 x 3 integer, the origin of the matrix
%               default is round(mat_123/2)
%
% Output:
%    The image file is stored as designated by filename.
%    
% written by Jinsoo Uh (jinsoo.uh@utsouthwestern.edu)
% 2007-03-01

function f = write_ANALYZE(varargin)

ImageArray = varargin{1};
filename = varargin{2};
mat_123 = varargin{3};
mat_123 = reshape(mat_123,1,3);

if nargin > 3 
   resolution = varargin{4};
   resolution = reshape(resolution,1,3);
else 
   resolution = [1 1 1];
end;

if nargin > 4
   scale_factor = varargin{5};
else
   scale_factor = 1;
end;

if nargin > 5
   data_type = varargin{6};
else
   data_type = 4;
end;

if nargin > 6
   no_header = varargin{7};
else
   no_header = 0;
end;
   
if nargin > 7
   mat_origin = varargin{8};
else
   mat_origin = round(mat_123/2);
end;

ImageArray = ImageArray * scale_factor;

[fname_path fname_body fname_ext] = fileparts(filename);
   fid_scan = fopen(fullfile(fname_path, [fname_body '.img']),'w'); 

      if data_type == 4
         fwrite(fid_scan,reshape(ImageArray,mat_123(1)*mat_123(2)*mat_123(3),1),'int16');
      elseif data_type == 2
         fwrite(fid_scan,reshape(ImageArray,mat_123(1)*mat_123(2)*mat_123(3),1),'uint8');
      elseif data_type == 16
         fwrite(fid_scan,reshape(ImageArray,mat_123(1)*mat_123(2)*mat_123(3),1),'float');
      elseif data_type == 64
         fwrite(fid_scan,reshape(ImageArray,mat_123(1)*mat_123(2)*mat_123(3),1),'double');
      end;          
    fclose(fid_scan);

    if no_header == 0
    P = fullfile(fname_path, [fname_body '.hdr']);

    DIM = mat_123;
    VOX = resolution;
    SCALE = 1;
    TYPE = data_type;
    OFFSET = 0;
    ORIGIN = mat_origin;
    DESCRIP = '';
    write_analyze75_header(P,DIM,VOX,SCALE,TYPE,OFFSET,ORIGIN,DESCRIP);
    end;
    

