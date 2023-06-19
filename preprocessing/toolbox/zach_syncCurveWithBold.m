function boldSyncedCurve = zach_syncCurveWithBold(envelopeCo2FilePath,globalDelay)

% writes a new file with altered timestamps. Global delay is added to all
% of the timestamps, ideally making t=0 the beginning of the curve synced
% with the beginning of the global BOLD signal

r1 = textread(envelopeCo2FilePath,'%f','delimiter',',');
timecourse = reshape(r1,2,length(r1)/2)';

% recall: negative delays should shift the etco2 curve to the left,
% positive delays shuld shift the etco2 curve to the right
timecourse(:,1) = timecourse(:,1)+globalDelay;

boldSyncedcurve = timecourse;

boldSyncedCurve = [envelopeCo2FilePath(1:end-4) '_boldSynced.txt'];

dlmwrite(boldSyncedCurve, boldSyncedcurve);

