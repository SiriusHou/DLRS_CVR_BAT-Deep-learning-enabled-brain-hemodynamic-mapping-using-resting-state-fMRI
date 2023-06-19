function etco2timecourse = zach_extractEtCO2(etco2FilePath)

% takes a path to a text file containing etco2 data, and returns a
% two-column vector with the timestamps in the first column and the co2
% values in the second column

%%
fid = fopen(etco2FilePath);
r1 = textscan(fid,'%f,%f','delimiter',',');
fclose(fid);
etco2timecourse = [r1{1},r1{2}];