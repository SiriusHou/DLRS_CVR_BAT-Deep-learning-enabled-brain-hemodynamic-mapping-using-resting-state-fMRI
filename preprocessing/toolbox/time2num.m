function Num = time2num(TimeString)
% Convert time format 'HH:MM:SS' to number (Unit: seconds)

t_vec = datevec(TimeString);
Num = t_vec(4)*3600 + t_vec(5)*60 + t_vec(6);
end