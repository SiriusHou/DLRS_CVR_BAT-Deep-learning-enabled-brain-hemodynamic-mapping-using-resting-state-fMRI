clc;
clear all;
close all;


%% Navigate to code folder and data folder
[code_directory, ~, ~] = fileparts(mfilename('fullpath'));
addpath(code_directory);

[parent_directory, ~, ~] = fileparts(code_directory);
data_directory = [parent_directory, filesep, 'data'];

subfolder = dir(data_directory);
subfolder(ismember({subfolder.name}, {'.', '..'})) =[];

%% Preprocessing pipeline
for ii=1:length(subfolder)
    RS_preprocessing(data_directory, subfolder(ii).name, code_directory);
    RS_BAT(data_directory, subfolder(ii).name, code_directory);
    RS_CVR_corrMap(data_directory, subfolder(ii).name, code_directory);
end
