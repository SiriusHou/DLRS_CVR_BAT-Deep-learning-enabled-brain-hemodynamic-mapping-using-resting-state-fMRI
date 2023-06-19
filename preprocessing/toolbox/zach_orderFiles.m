function orderedResult = zach_orderFiles(outOrdered)
%%
% takes a character array where each line is a file. This program must
% order the character array in the following manner:

% filepath\myfile_1zb.txt
% filepath\myfile_2zb.txt
% filepath\myfile_3zb.txt
% filepath\myfile_4zb.txt
% filepath\myfile_5zb.txt
% ...
% filepath\myfile_19zb.txt
% filepath\myfile_20zb.txt
% filepath\myfile_21zb.txt
% etc.

%%

% convert char array to cells containing each line
line = '';
lines = cell(size(outOrdered,1),1);
names = cell(size(outOrdered,1),1);
for i = 1:size(outOrdered,1)
    lines{i} = outOrdered(i,:);
    [~,name] = fileparts(outOrdered(i,:));
    names{i} = name;
end

data = [names,lines];

orderedResult = lines;
