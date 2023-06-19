function ASL_downSampleMNI(fname,outfile)
% Resample the MNI space matrix [181 217 181] resolution [1 1 1] to matrix
% [91 109 91] resolution [2 2 2] origin [46 64 37]
% example: 
%   ASL_downSampleMNI('C:\Results\test_2d_asl_aCBF_thre_MNI.img')

% [pathstr,name,ext] = fileparts(fname);

old_img = loadimage(fname);
new_img = affine(old_img, [1 0 0 0;0 1 0 0;0 0 1 0;0 0 0 1], [2 2 2], 0);

write_ANALYZE(new_img,outfile,[91 109 91],[2 2 2],1,16,0,[46 64 37]);


function [new_img, new_M] = affine(old_img, old_M, new_elem_size, verbose, bg, method)
%  Using 2D or 3D affine matrix to rotate, translate, scale, reflect and
%  shear a 2D image or 3D volume. 2D image is represented by a 2D matrix,
%  3D volume is represented by a 3D matrix, and data type can be real 
%  integer or floating-point.
%
%  You may notice that MATLAB has a function called 'imtransform.m' for
%  2D spatial transformation. However, keep in mind that 'imtransform.m'
%  assumes y for the 1st dimension, and x for the 2nd dimension. They are
%  equivalent otherwise.
%
%  In addition, if you adjust the 'new_elem_size' parameter, this 'affine.m'
%  is equivalent to 'interp2.m' for 2D image, and equivalent to 'interp3.m'
%  for 3D volume.
%
%  Usage: [new_img new_M] = ...
%	affine(old_img, old_M, [new_elem_size], [verbose], [bg], [method]);
%
%  old_img  -	original 2D image or 3D volume. We assume x for the 1st
%		dimension, y for the 2nd dimension, and z for the 3rd
%		dimension.
%
%  old_M  -	a 3x3 2D affine matrix for 2D image, or a 4x4 3D affine
%		matrix for 3D volume. We assume x for the 1st dimension,
%		y for the 2nd dimension, and z for the 3rd dimension.
%
%  new_elem_size (optional)  -  size of voxel along x y z direction for 
%		a transformed 3D volume, or size of pixel along x y for
%		a transformed 2D image. We assume x for the 1st dimension
%		y for the 2nd dimension, and z for the 3rd dimension.
%		'new_elem_size' is 1 if it is default or empty.
%
%		You can increase its value to decrease the resampling rate,
%		and make the 2D image or 3D volume more coarse. It works
%		just like 'interp3'.
%
%  verbose (optional) - 1, 0
%		1:  show transforming progress in percentage
%		2:  progress will not be displayed
%		'verbose' is 1 if it is default or empty.
%
%  bg (optional)  -	background voxel intensity in any extra corner that
%		is caused by the interpolation. 0 in most cases. If it is
%		default or empty, 'bg' will be the average of two corner
%		voxel intensities in original data.
%
%  method (optional)  -	1, 2, or 3
%		1:  for Trilinear interpolation
%		2:  for Nearest Neighbor interpolation
%		3:  for Fischer's Bresenham interpolation
%		'method' is 1 if it is default or empty.
%
%  new_img  -	transformed 2D image or 3D volume
%
%  new_M  -	transformed affine matrix
%
%  Example 1 (3D rotation):
%	load mri.mat;   old_img = double(squeeze(D));
%	old_M = [0.88 0.5 3 -90; -0.5 0.88 3 -126; 0 0 2 -72; 0 0 0 1];
%	new_img = affine(old_img, old_M, 2);
%	[x y z] = meshgrid(1:128,1:128,1:27);
%	sz = size(new_img);
%	[x1 y1 z1] = meshgrid(1:sz(2),1:sz(1),1:sz(3));
%	figure; slice(x, y, z, old_img, 64, 64, 13.5);
%	shading flat; colormap(map); view(-66, 66);
%	figure; slice(x1, y1, z1, new_img, sz(1)/2, sz(2)/2, sz(3)/2);
%	shading flat; colormap(map); view(-66, 66);
%
%  Example 2 (2D interpolation):
%	load mri.mat;   old_img=D(:,:,1,13)';
%	old_M = [1 0 0; 0 1 0; 0 0 1];
%	new_img = affine(old_img, old_M, [.2 .4]);
%	figure; image(old_img); colormap(map);
%	figure; image(new_img); colormap(map);
%
%  This program is inspired by:
%  SPM5 Software from Wellcome Trust Centre for Neuroimaging
%	http://www.fil.ion.ucl.ac.uk/spm/software
%  Fischer, J., A. del Rio (2004). A Fast Method for Applying Rigid
%	Transformations to Volume Data, WSCG2004 Conference.
%	http://wscg.zcu.cz/wscg2004/Papers_2004_Short/M19.pdf
%  
%  - Jimmy Shen (jimmy@rotman-baycrest.on.ca)

   if ~exist('old_img','var') | ~exist('old_M','var')
      error('Usage: [new_img new_M] = affine(old_img, old_M, [new_elem_size], [verbose], [bg], [method]);');
   end

   if ndims(old_img) == 3
      if ~isequal(size(old_M),[4 4])
         error('old_M should be a 4x4 affine matrix for 3D volume.');
      end
   elseif ndims(old_img) == 2
      if ~isequal(size(old_M),[3 3])
         error('old_M should be a 3x3 affine matrix for 2D image.');
      end
   else
      error('old_img should be either 2D image or 3D volume.');
   end

   if ~exist('new_elem_size','var') | isempty(new_elem_size)
      new_elem_size = [1 1 1];
   elseif length(new_elem_size) < 2
      new_elem_size = new_elem_size(1)*ones(1,3);
   elseif length(new_elem_size) < 3
      new_elem_size = [new_elem_size(:); 1]';
   end

   if ~exist('method','var') | isempty(method)
      method = 1;
   elseif ~exist('bresenham_line3d.m','file') & method == 3
      error([char(10) char(10) 'Please download 3D Bresenham''s line generation program from:' char(10) char(10) 'http://www.mathworks.com/matlabcentral/fileexchange/loadFile.do?objectId=21057' char(10) char(10) 'to test Fischer''s Bresenham interpolation method.' char(10) char(10)]);
   end

   %  Make compatible to MATLAB earlier than version 7 (R14), which
   %  can only perform arithmetic on double data type
   %
   old_img = double(old_img);
   old_dim = size(old_img);

   if ~exist('bg','var') | isempty(bg)
      bg = mean([old_img(1) old_img(end)]);
   end

   if ~exist('verbose','var') | isempty(verbose)
      verbose = 1;
   end

   if ndims(old_img) == 2
      old_dim(3) = 1;
      old_M = old_M(:, [1 2 3 3]);
      old_M = old_M([1 2 3 3], :);
      old_M(3,:) = [0 0 1 0];
      old_M(:,3) = [0 0 1 0]';
   end

   %  Vertices of img in voxel
   %
   XYZvox = [	1		1		1
		1		1		old_dim(3)
		1		old_dim(2)	1
		1		old_dim(2)	old_dim(3)
		old_dim(1)	1		1
		old_dim(1)	1		old_dim(3)
		old_dim(1)	old_dim(2)	1
		old_dim(1)	old_dim(2)	old_dim(3)   ]';

   old_R = old_M(1:3,1:3);
   old_T = old_M(1:3,4);

   %  Vertices of img in millimeter
   %
   XYZmm = old_R*(XYZvox-1) + repmat(old_T, [1, 8]);

   %  Make scale of new_M according to new_elem_size
   %
   new_M = diag([new_elem_size 1]);

   %  Make translation so minimum vertex is moved to [1,1,1]
   %
   new_M(1:3,4) = round( min(XYZmm,[],2) );

   %  New dimensions will be the maximum vertices in XYZ direction (dim_vox)
   %  i.e. compute   dim_vox   via   dim_mm = R*(dim_vox-1)+T
   %  where, dim_mm = round(max(XYZmm,[],2));
   %
   new_dim = ceil(new_M(1:3,1:3) \ ( round(max(XYZmm,[],2))-new_M(1:3,4) )+1)';

   %  Initialize new_img with new_dim
   %
   new_img = zeros(new_dim(1:3));

   %  Mask out any changes from Z axis of transformed volume, since we
   %  will traverse it voxel by voxel below. We will only apply unit
   %  increment of mask_Z(3,4) to simulate the cursor movement
   %
   %  i.e. we will use   mask_Z * new_XYZvox   to replace   new_XYZvox
   %
   mask_Z = diag(ones(1,4));
   mask_Z(3,3) = 0;

   %  It will be easier to do the interpolation if we invert the process
   %  by not traversing the original volume. Instead, we traverse the
   %  transformed volume, and backproject each voxel in the transformed 
   %  volume back into the original volume. If the backprojected voxel
   %  in original volume is within its boundary, the intensity of that
   %  voxel can be used by the cursor location in the transformed volume.
   %
   %  First, we traverse along Z axis of transformed volume voxel by voxel
   %
   for z = 1:new_dim(3)

      if verbose & ~mod(z,10)
         fprintf('%.2f percent is done.\n', 100*z/new_dim(3));
      end

      %  We need to find out the mapping from voxel in the transformed
      %  volume (new_XYZvox) to voxel in the original volume (old_XYZvox)
      %
      %  The following equation works, because they all equal to XYZmm:
      %  new_R*(new_XYZvox-1) + new_T  ==  old_R*(old_XYZvox-1) + old_T
      %
      %  We can use modified new_M1 & old_M1 to substitute new_M & old_M
      %      new_M1 * new_XYZvox       ==       old_M1 * old_XYZvox
      %
      %  where: M1 = M;   M1(:,4) = M(:,4) - sum(M(:,1:3),2);
      %  and:             M(:,4) == [T; 1] == sum(M1,2)
      %
      %  Therefore:   old_XYZvox = old_M1 \ new_M1 * new_XYZvox;
      %
      %  Since we are traverse Z axis, and   new_XYZvox   is replaced
      %  by   mask_Z * new_XYZvox, the above formula can be rewritten
      %  as:    old_XYZvox = old_M1 \ new_M1 * mask_Z * new_XYZvox;
      %
      %  i.e. we find the mapping from new_XYZvox to old_XYZvox:
      %  M = old_M1 \ new_M1 * mask_Z;
      %
      %  First, compute modified old_M1 & new_M1
      %
      old_M1 = old_M;   old_M1(:,4) = old_M(:,4) - sum(old_M(:,1:3),2);
      new_M1 = new_M;   new_M1(:,4) = new_M(:,4) - sum(new_M(:,1:3),2);

      %  Then, apply unit increment of mask_Z(3,4) to simulate the
      %  cursor movement
      %
      mask_Z(3,4) = z;

      %  Here is the mapping from new_XYZvox to old_XYZvox
      %
      M = old_M1 \ new_M1 * mask_Z;

      switch method
      case 1
         new_img(:,:,z) = trilinear(old_img, new_dim, old_dim, M, bg);
      case 2
         new_img(:,:,z) = nearest_neighbor(old_img, new_dim, old_dim, M, bg);
      case 3
         new_img(:,:,z) = bresenham(old_img, new_dim, old_dim, M, bg);
      end

   end;			% for z

   if ndims(old_img) == 2
      new_M(3,:) = [];
      new_M(:,3) = [];
   end

   return;					% affine


%--------------------------------------------------------------------
function img_slice = trilinear(img, dim1, dim2, M, bg)

   img_slice = zeros(dim1(1:2));
   TINY = 5e-2;					% tolerance

   %  Dimension of transformed 3D volume
   %
   xdim1 = dim1(1);
   ydim1 = dim1(2);

   %  Dimension of original 3D volume
   %
   xdim2 = dim2(1);
   ydim2 = dim2(2);
   zdim2 = dim2(3);

   %  initialize new_Y accumulation
   %
   Y2X = 0;
   Y2Y = 0;
   Y2Z = 0;

   for y = 1:ydim1

      %  increment of new_Y accumulation
      %
      Y2X = Y2X + M(1,2);		% new_Y to old_X
      Y2Y = Y2Y + M(2,2);		% new_Y to old_Y
      Y2Z = Y2Z + M(3,2);		% new_Y to old_Z

      %  backproject new_Y accumulation and translation to old_XYZ
      %
      old_X = Y2X + M(1,4);
      old_Y = Y2Y + M(2,4);
      old_Z = Y2Z + M(3,4);

      for x = 1:xdim1

         %  accumulate the increment of new_X, and apply it
         %  to the backprojected old_XYZ
         %
         old_X = M(1,1) + old_X  ;
         old_Y = M(2,1) + old_Y  ;
         old_Z = M(3,1) + old_Z  ;

         %  within boundary of original image
         %
         if (	old_X > 1-TINY & old_X < xdim2+TINY & ...
		old_Y > 1-TINY & old_Y < ydim2+TINY & ...
		old_Z > 1-TINY & old_Z < zdim2+TINY	)

            %  Calculate distance of old_XYZ to its neighbors for
            %  weighted intensity average
            %
            dx = old_X - floor(old_X);
            dy = old_Y - floor(old_Y);
            dz = old_Z - floor(old_Z);

            x000 = floor(old_X);
            x100 = x000 + 1;

            if floor(old_X) < 1
               x000 = 1;
               x100 = x000;
            elseif floor(old_X) > xdim2-1
               x000 = xdim2;
               x100 = x000;
            end

            x010 = x000;
            x001 = x000;
            x011 = x000;

            x110 = x100;
            x101 = x100;
            x111 = x100;

            y000 = floor(old_Y);
            y010 = y000 + 1;

            if floor(old_Y) < 1
               y000 = 1;
               y100 = y000;
            elseif floor(old_Y) > ydim2-1
               y000 = ydim2;
               y010 = y000;
            end

            y100 = y000;
            y001 = y000;
            y101 = y000;

            y110 = y010;
            y011 = y010;
            y111 = y010;

            z000 = floor(old_Z);
            z001 = z000 + 1;

            if floor(old_Z) < 1
               z000 = 1;
               z001 = z000;
            elseif floor(old_Z) > zdim2-1
               z000 = zdim2;
               z001 = z000;
            end

            z100 = z000;
            z010 = z000;
            z110 = z000;

            z101 = z001;
            z011 = z001;
            z111 = z001;

            x010 = x000;
            x001 = x000;
            x011 = x000;

            x110 = x100;
            x101 = x100;
            x111 = x100;

            v000 = double(img(x000, y000, z000));
            v010 = double(img(x010, y010, z010));
            v001 = double(img(x001, y001, z001));
            v011 = double(img(x011, y011, z011));

            v100 = double(img(x100, y100, z100));
            v110 = double(img(x110, y110, z110));
            v101 = double(img(x101, y101, z101));
            v111 = double(img(x111, y111, z111));

            img_slice(x,y) = v000*(1-dx)*(1-dy)*(1-dz) + ...
               v010*(1-dx)*dy*(1-dz) + ...
               v001*(1-dx)*(1-dy)*dz + ...
               v011*(1-dx)*dy*dz + ...
               v100*dx*(1-dy)*(1-dz) + ...
               v110*dx*dy*(1-dz) + ...
               v101*dx*(1-dy)*dz + ...
               v111*dx*dy*dz;

         else
            img_slice(x,y) = bg;

         end	% if boundary

      end	% for x
   end		% for y

   return;					% trilinear


%--------------------------------------------------------------------
function img_slice = nearest_neighbor(img, dim1, dim2, M, bg)

   img_slice = zeros(dim1(1:2));

   %  Dimension of transformed 3D volume
   %
   xdim1 = dim1(1);
   ydim1 = dim1(2);

   %  Dimension of original 3D volume
   %
   xdim2 = dim2(1);
   ydim2 = dim2(2);
   zdim2 = dim2(3);

   %  initialize new_Y accumulation
   %
   Y2X = 0;
   Y2Y = 0;
   Y2Z = 0;

   for y = 1:ydim1

      %  increment of new_Y accumulation
      %
      Y2X = Y2X + M(1,2);		% new_Y to old_X
      Y2Y = Y2Y + M(2,2);		% new_Y to old_Y
      Y2Z = Y2Z + M(3,2);		% new_Y to old_Z

      %  backproject new_Y accumulation and translation to old_XYZ
      %
      old_X = Y2X + M(1,4);
      old_Y = Y2Y + M(2,4);
      old_Z = Y2Z + M(3,4);

      for x = 1:xdim1

         %  accumulate the increment of new_X and apply it
         %  to the backprojected old_XYZ
         %
         old_X = M(1,1) + old_X  ;
         old_Y = M(2,1) + old_Y  ;
         old_Z = M(3,1) + old_Z  ;

         xi = round(old_X);
         yi = round(old_Y);
         zi = round(old_Z);

         %  within boundary of original image
         %
         if (	xi >= 1 & xi <= xdim2 & ...
		yi >= 1 & yi <= ydim2 & ...
		zi >= 1 & zi <= zdim2	)

            img_slice(x,y) = img(xi,yi,zi);

         else
            img_slice(x,y) = bg;

         end	% if boundary

      end	% for x
   end		% for y

   return;					% nearest_neighbor


%--------------------------------------------------------------------
function img_slice = bresenham(img, dim1, dim2, M, bg)

   img_slice = zeros(dim1(1:2));

   %  Dimension of transformed 3D volume
   %
   xdim1 = dim1(1);
   ydim1 = dim1(2);

   %  Dimension of original 3D volume
   %
   xdim2 = dim2(1);
   ydim2 = dim2(2);
   zdim2 = dim2(3);

   for y = 1:ydim1

      start_old_XYZ = round(M*[0     y 0 1]');
      end_old_XYZ   = round(M*[xdim1 y 0 1]');

      [X Y Z] = bresenham_line3d(start_old_XYZ, end_old_XYZ);

      %  line error correction
      %
%      del = end_old_XYZ - start_old_XYZ;
 %     del_dom = max(del);
  %    idx_dom = find(del==del_dom);
   %   idx_dom = idx_dom(1);
    %  idx_other = [1 2 3];
     % idx_other(idx_dom) = [];
      %del_x1 = del(idx_other(1));
%      del_x2 = del(idx_other(2));
 %     line_slope = sqrt((del_x1/del_dom)^2 + (del_x2/del_dom)^2 + 1);
  %    line_error = line_slope - 1;
% line error correction removed because it is too slow

      for x = 1:xdim1

         %  rescale ratio
         %
         i = round(x * length(X) / xdim1);

         if i < 1
            i = 1;
         elseif i > length(X)
            i = length(X);
         end

         xi = X(i);
         yi = Y(i);
         zi = Z(i);

         %  within boundary of the old XYZ space
         %
         if (	xi >= 1 & xi <= xdim2 & ...
		yi >= 1 & yi <= ydim2 & ...
		zi >= 1 & zi <= zdim2	)

            img_slice(x,y) = img(xi,yi,zi);

%            if line_error > 1
 %              x = x + 1;

%               if x <= xdim1
 %                 img_slice(x,y) = img(xi,yi,zi);
  %                line_error = line_slope - 1;
   %            end
    %        end		% if line_error
% line error correction removed because it is too slow

         else
            img_slice(x,y) = bg;

         end	% if boundary

      end	% for x
   end		% for y

   return;					% bresenham



function [f varargout] = loadimage(varargin)
% This function read a volume of image without needing to enter 
% the data type and the matrix size. Those information is read
% from the header file. 
%   Input:  filename ~ the file name e.g.) 'mpr.img'
%           l_scale  ~ if 1 multiply the scale in .hdr file
%                      default value is 0         
%      You may use wildcards (i.e., * ) in the file name for the image files
%        of the same structure.
%      Then, the volumes of the images will be loaded to 4D matrix
%      in the sequence of the file name order.
%
%   Output : the 3D or 4D array of the image
%
     filename = varargin{1};
     if nargin > 1
        l_scale = varargin{2};
     else
        l_scale = 0;
     end;
        
     filelist = dir(filename);
     file_path=fileparts(filename);

     [DIM,VOX,SCALE,TYPE,OFFSET,ORIGIN,DESCRIP] = read_analyze75_header(fullfile(file_path,filelist(1,1).name));
     array = zeros(DIM(1),DIM(2),DIM(3),length(filelist));
     

for fileindex = 1:length(filelist)

     [DIM,VOX,SCALE,TYPE,OFFSET,ORIGIN,DESCRIP] = read_analyze75_header(fullfile(file_path,filelist(fileindex,1).name));    

     fid = fopen(fullfile(file_path,filelist(fileindex,1).name),'r');
     switch TYPE
     case 2,
         dummy = fread(fid, DIM(1)*DIM(2)*DIM(3),'uint8');
     case 4,
         dummy = fread(fid, DIM(1)*DIM(2)*DIM(3),'int16');
     case 16,
         dummy = fread(fid, DIM(1)*DIM(2)*DIM(3),'float');
     case 64,
         dummy = fread(fid, DIM(1)*DIM(2)*DIM(3),'double');
     end;
     fclose(fid);
     
     if l_scale == 1
        array(:,:,:,fileindex) = reshape(dummy, DIM)*SCALE;
     else
        array(:,:,:,fileindex) = reshape(dummy, DIM);
     end;  
   
end;
    
     
     f = array;
     varargout{1} = DIM;
     varargout{2} = VOX;
     varargout{3} = SCALE;     
     varargout{4} = TYPE;          
     varargout{5} = OFFSET;          
     varargout{6} = ORIGIN;          
     
function [DIM,VOX,SCALE,TYPE,OFFSET,ORIGIN,DESCRIP] = read_analyze75_header(P)
% reads a header
% FORMAT [DIM VOX SCALE TYPE OFFSET ORIGIN DESCRIP] = spm_hread(P);
%
% P       - filename 	     (e.g spm or spm.img)
% DIM     - image size       [i j k [l]] (voxels)
% VOX     - voxel size       [x y z [t]] (mm [secs])
% SCALE   - scale factor
% TYPE    - datatype (integer - see spm_type)
% OFFSET  - offset (bytes)
% ORIGIN  - origin [i j k]
% DESCRIP - description string
%___________________________________________________________________________
%
% spm_hread reads variables into working memory from a SPM/ANALYZE
% compatible header file.  If the header does not exist global defaults
% are used.  The 'originator' field of the ANALYZE format has been
% changed to ORIGIN in the SPM version of the header.  funused1 of the
% ANALYZE format is used for SCALE
%
% see also dbh.h (ANALYZE) spm_hwrite.m and spm_type.m
%
%__________________________________________________________________________
% @(#)spm_hread.m	2.7 99/10/29

% ensure correct suffix {.hdr}
%---------------------------------------------------------------------------
P     = deblank(P);
q     = length(P);
if q>=4 & P(q - 3) == '.'; P = P(1:(q - 4)); end
P     = [P '.hdr'];

% open header file
%---------------------------------------------------------------------------
fid   = fopen(P,'r','native');

if (fid > 0)

% read (struct) header_key
%---------------------------------------------------------------------------
fseek(fid,0,'bof');

otherendian = 0;
sizeof_hdr 	= fread(fid,1,'int32');
if sizeof_hdr==1543569408, % Appears to be other-endian
	% Re-open other-endian
	fclose(fid);
	if spm_platform('bigend'),
		fid = fopen(P,'r','ieee-le');
	else,
		fid = fopen(P,'r','ieee-be');
	end;
	fseek(fid,0,'bof');
	sizeof_hdr = fread(fid,1,'int32');
	otherendian = 1;
end;

data_type  	= mysetstr(fread(fid,10,'uchar'))';
db_name    	= mysetstr(fread(fid,18,'uchar'))';
extents    	= fread(fid,1,'int32');
session_error   = fread(fid,1,'int16');
regular    	= mysetstr(fread(fid,1,'uchar'))';
hkey_un0    	= mysetstr(fread(fid,1,'uchar'))';



% read (struct) image_dimension
%---------------------------------------------------------------------------
fseek(fid,40,'bof');

dim    		= fread(fid,8,'int16');
vox_units    	= mysetstr(fread(fid,4,'uchar'))';
cal_units    	= mysetstr(fread(fid,8,'uchar'))';
unused1		= fread(fid,1,'int16');
datatype	= fread(fid,1,'int16');
bitpix		= fread(fid,1,'int16');
dim_un0		= fread(fid,1,'int16');
pixdim		= fread(fid,8,'float');
vox_offset	= fread(fid,1,'float');
funused1	= fread(fid,1,'float');
funused2	= fread(fid,1,'float');
funused3	= fread(fid,1,'float');
cal_max		= fread(fid,1,'float');
cal_min		= fread(fid,1,'float');
compressed	= fread(fid,1,'int32');
verified	= fread(fid,1,'int32');
glmax		= fread(fid,1,'int32');
glmin		= fread(fid,1,'int32');

% read (struct) data_history
%---------------------------------------------------------------------------
fseek(fid,148,'bof');

descrip		= mysetstr(fread(fid,80,'uchar'))';
aux_file	= mysetstr(fread(fid,24,'uchar'))';
orient		= fread(fid,1,'uchar');
origin		= fread(fid,5,'int16');
generated	= mysetstr(fread(fid,10,'uchar'))';
scannum		= mysetstr(fread(fid,10,'uchar'))';
patient_id	= mysetstr(fread(fid,10,'uchar'))';
exp_date	= mysetstr(fread(fid,10,'uchar'))';
exp_time	= mysetstr(fread(fid,10,'uchar'))';
hist_un0	= mysetstr(fread(fid,3,'uchar'))';
views		= fread(fid,1,'int32');
vols_added	= fread(fid,1,'int32');
start_field	= fread(fid,1,'int32');
field_skip	= fread(fid,1,'int32');
omax		= fread(fid,1,'int32');
omin		= fread(fid,1,'int32');
smax		= fread(fid,1,'int32');
smin		= fread(fid,1,'int32');

fclose(fid);

if isempty(smin)
	error(['There is a problem with the header file ' P '.']);
end

% convert to SPM global variables
%---------------------------------------------------------------------------
DIM    	  	= dim(2:4)';
VOX       	= pixdim(2:4)';
SCALE     	= funused1;
SCALE    	= ~SCALE + SCALE;
TYPE     	= datatype;
if otherendian == 1 & datatype ~= 2,
	TYPE = TYPE*256;
end;
OFFSET	  	= vox_offset;
ORIGIN    	= origin(1:3)';
DESCRIP   	= descrip(1:max(find(descrip)));

else
	global DIM VOX SCALE TYPE OFFSET ORIGIN
	DESCRIP = ['defaults'];
end
return;
%_______________________________________________________________________

function out = mysetstr(in)
tmp = find(in == 0);
tmp = min([min(tmp) length(in)]);
out = setstr(in(1:tmp));
return;


function f = write_ANALYZE(varargin)
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
    

function f = write_analyze75_header(P,DIM,VOX,SCALE,TYPE,OFFSET,ORIGIN,DESCRIP)
% writes an ANALYZE header 
% modified from spm_hwrite.m
%
% P       - filename 	     (e.g 'spm' or 'spm.img')
% DIM     - image size       [i j k [l]] (voxels)
% VOX     - voxel size       [x y z [t]] (mm [sec])
% SCALE   - scale factor
% TYPE    - datatype (integer - see spm_type)
% OFFSET  - offset (bytes)
% ORIGIN  - [i j k] of origin  (default = [0 0 0])
% DESCRIP - description string (default = 'spm compatible')
%

P               = P(P ~= ' ');
q    		= length(P);
if q>=4 & P(q - 3) == '.', P = P(1:(q - 4)); end;
P     		= [P '.hdr'];

fid             = fopen(P,'w','native');
%fid             = fopen(P,'w','ieee-le');
%fid             = fopen(P,'w','ieee-be');

%---------------------------------------------------------------------------
data_type 	= ['dsr      ' 0];

P     		= [P '                  '];
db_name		= [P(1:17) 0];

% set header variables
%---------------------------------------------------------------------------
DIM		= DIM(:)'; if size(DIM,2) < 4; DIM = [DIM 1]; end
VOX		= VOX(:)'; if size(VOX,2) < 4; VOX = [VOX 0]; end
dim		= [4 DIM(1:4) 0 0 0];	
pixdim		= [0 VOX(1:4) 0 0 0];
vox_offset      = OFFSET;
funused1	= SCALE;
glmax		= 1;
glmin		= 0;
bitpix 		= 0;
descrip         = zeros(1,80);
aux_file        = ['none                   ' 0];
origin          = [0 0 0 0 0];

%---------------------------------------------------------------------------
if TYPE == 1;   bitpix = 1;  glmax = 1;        glmin = 0;	end
if TYPE == 2;   bitpix = 8;  glmax = 255;      glmin = 0;	end
if TYPE == 4;   bitpix = 16; glmax = 32767;    glmin = 0;  	end
if TYPE == 8;   bitpix = 32; glmax = (2^31-1); glmin = 0;	end
if TYPE == 16;  bitpix = 32; glmax = 1;        glmin = 0;	end
if TYPE == 64;  bitpix = 64; glmax = 1;        glmin = 0;	end

%---------------------------------------------------------------------------
if nargin >= 7; origin = [ORIGIN(:)' 0 0];  end
if nargin <  8; DESCRIP = 'spm compatible'; end

d          	= 1:min([length(DESCRIP) 79]);
descrip(d) 	= DESCRIP(d);

fseek(fid,0,'bof');

% write (struct) header_key
%---------------------------------------------------------------------------
fwrite(fid,348,		'int32');
fwrite(fid,data_type,	'char' );
fwrite(fid,db_name,	'char' );
fwrite(fid,0,		'int32');
fwrite(fid,0,		'int16');
fwrite(fid,'r',		'char' );
fwrite(fid,'0',		'char' );

% write (struct) image_dimension
%---------------------------------------------------------------------------
fseek(fid,40,'bof');

fwrite(fid,dim,		'int16');
fwrite(fid,'mm',	'char' );
fwrite(fid,0,		'char' );
fwrite(fid,0,		'char' );

fwrite(fid,zeros(1,8),	'char' );
fwrite(fid,0,		'int16');
fwrite(fid,TYPE,	'int16');
fwrite(fid,bitpix,	'int16');
fwrite(fid,0,		'int16');
fwrite(fid,pixdim,	'float');
fwrite(fid,vox_offset,	'float');
fwrite(fid,funused1,	'float');
fwrite(fid,0,		'float');
fwrite(fid,0,		'float');
fwrite(fid,0,		'float');
fwrite(fid,0,		'float');
fwrite(fid,0,		'int32');
fwrite(fid,0,		'int32');
fwrite(fid,glmax,	'int32');
fwrite(fid,glmin,	'int32');

% write (struct) image_dimension
%---------------------------------------------------------------------------
fwrite(fid,descrip,	'char');
fwrite(fid,aux_file,    'char');
fwrite(fid,0,           'char');
fwrite(fid,origin,      'int16');

if fwrite(fid,zeros(1,85), 'char')~=85
	fclose(fid);
end

fclose(fid);


