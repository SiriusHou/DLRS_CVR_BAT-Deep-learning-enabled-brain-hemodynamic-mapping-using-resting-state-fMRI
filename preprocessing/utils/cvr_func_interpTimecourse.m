function etco2curve = cvr_func_interpTimecourse(etco2timecourse,shift,dynnum,timestep)
% interpolate timecourse ([time_col, value_col]) according to shift,
% dynnum, and timesteps.
% if no value, fill with nearest available value in timecourse


maxtime = etco2timecourse(:,1);
maxco2  = etco2timecourse(:,2);

lastTime    = max(maxtime);
firstTime   = min(maxtime);
z           = -shift + (0:(dynnum-1))*timestep;
    
% find if any time values are outside of timecourse range
ind1 = find(z <  firstTime);
ind2 = find((z >= firstTime) & (z <= lastTime));
ind3 = find(z >  lastTime);

% set all time values < starting time to be the first value
if(~isempty(ind1))
    etco2curve(ind1) = maxco2(1);
end;
% set all time values > ending time to be the last value.
if(~isempty(ind3))
    etco2curve(ind3) = maxco2(end);
end;

etco2curve(ind2) = interp1(maxtime,maxco2,z(ind2),'linear');
etco2curve = etco2curve';


