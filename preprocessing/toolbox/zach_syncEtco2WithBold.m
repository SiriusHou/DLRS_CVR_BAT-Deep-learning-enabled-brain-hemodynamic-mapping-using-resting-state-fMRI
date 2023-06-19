function boldSyncedEtco2Path = zach_syncEtco2WithBold(envelopeCo2FilePath,globalDelay)

% writes a new file with altered timestamps. Global delay is added to all
% of the timestamps, ideally making t=0 the beginning of the EtCO2 synced
% with the beginning of the global BOLD signal

r1 = textread(envelopeCo2FilePath,'%f','delimiter',',');
etco2timecourse = reshape(r1,2,length(r1)/2)';

% recall: negative delays should shift the etco2 curve to the left,
% positive delays shuld shift the etco2 curve to the right
etco2timecourse(:,1) = etco2timecourse(:,1)+globalDelay;

boldSyncedEtco2 = etco2timecourse;

boldSyncedEtco2Path = [envelopeCo2FilePath(1:end-4) '_boldSynced.txt'];

dlmwrite(boldSyncedEtco2Path, boldSyncedEtco2);

