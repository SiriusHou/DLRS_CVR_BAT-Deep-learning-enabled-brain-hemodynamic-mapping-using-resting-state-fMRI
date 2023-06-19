function eto2timecourse = zach_extractEtO2Data(eto2FilePath)

% takes a path to a text file containing eto2 data, and returns a
% two-column vector with the timestamps in the first column and the co2
% values in the second column

fid = fopen(eto2FilePath);
data = textscan(fid,'%f %f','delimiter',{sprintf('\n'),','});
fclose(fid);
eto2timecourse = [data{1},data{2}];

end