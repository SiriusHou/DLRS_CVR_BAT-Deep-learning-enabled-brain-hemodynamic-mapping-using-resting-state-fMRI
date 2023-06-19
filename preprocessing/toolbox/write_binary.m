% ImageArray
% filename
% data_type

function f = write_binary(varargin)

ImageArray = varargin{1};
filename = varargin{2};
data_type = varargin{3};

fid = fopen(filename,'w'); 
fwrite(fid,ImageArray(:),data_type);
fclose(fid);
