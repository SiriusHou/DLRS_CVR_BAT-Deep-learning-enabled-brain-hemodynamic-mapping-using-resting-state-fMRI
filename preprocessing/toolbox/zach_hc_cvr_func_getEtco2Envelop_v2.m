function [txtname_etco2] = zach_hc_cvr_func_getEtco2Envelop_v2(physiodir,...
    txtname_co2_raw, sample_rate_co2, envelope_interp_rate, select)
% read etco2 recording and export envelop in text

% NO BOLD START TIME

disp('Finding EtCO2 Envelope...')
%% input overview
% physiodir           - directory of the co2 trace
% txtname_co2_raw     - name of the co2 trace text file
% time_bold_start     - time when co2 modulations begin
% bold_scandur        - duration of the experiment
% select              - control parameter in picking peak (0.1)
% tr                  - Repetition time (unit: sec)

% physiodir = 'C:\Users\Zach\OneDrive\College_Grad\Lu_Lab\peiyingYangStitchedPipeline\workspace';
% txtname_co2_raw = 'co2_recording.txt';
% select = 0.1;

%% inputs
fn_co2              = [physiodir filesep txtname_co2_raw];


%% deal with etco2 recording
% load text file
fid         = fopen(fn_co2);
data_co2    = textscan(fid,'%s %s','Delimiter',',\t','headerlines', 1); fclose(fid);
vec_co2     = str2num(char(data_co2{2}));

% 
timep_etco2 = ( 1:length(vec_co2) )' / sample_rate_co2;
value_etco2 = smooth(vec_co2,21);
curve_etco2 = [timep_etco2,value_etco2];

% % sort out the timing
% idx_1stpoint    = (timer_sync_mri + sample_delay_co2 - time_add_ahead) * sample_rate_co2;
% idx_mri_slot    = idx_1stpoint + (1 : (scan_duration * sample_rate_co2));
% curve_etco2     = [timep_etco2(idx_mri_slot), value_etco2(idx_mri_slot)];

% downsample average every 10 data points
span         = 10;
num_etco2    = size(curve_etco2,1);
num2         = floor(num_etco2/span) * span;
curve2_etco2 = [curve_etco2(floor(span/2):span:num2,1), ...
               mean(reshape(curve_etco2(1:num2,2),span,[]),1)'];


%% get the envelop of etco2 timecourse

% run yang's algorithm (finds peaks with scale-space algorithm)
timeline = curve2_etco2(:,1);
co2value = curve2_etco2(:,2);

[~,criterion]   = pickpeaks(co2value,select,0);
locs            = find(criterion > 0.01);

co2peaks_yang    = co2value(locs);
timepeaks_yang   = timeline(locs);

% run cuimei's algorithm
output = zach_cuimeiEtco2Alg(timeline,co2value);

timepeaks_cuimei = output(:,1);
co2peaks_cuimei = output(:,2);

% run 5 point median filter on both co2 curves
co2peaks_yang = medfilt1(co2peaks_yang,5);
co2peaks_cuimei = medfilt1(co2peaks_cuimei,5);

% find the latest start time and the earliest end time of the curves
if timepeaks_yang(1) > timepeaks_cuimei(1)
    latestStartTime = timepeaks_yang(1);
else
    latestStartTime = timepeaks_cuimei(1);
end
if timepeaks_yang(end) < timepeaks_cuimei(end)
    earliestEndTime = timepeaks_yang(end);
else
    earliestEndTime = timepeaks_cuimei(end);
end

% generate interploation timestamps (frequency in Hz specified by
% envelope_interp_rate)
interpTime = latestStartTime:1/envelope_interp_rate:earliestEndTime;

% interpolate curves
interpCo2PeaksYang = interp1(timepeaks_yang,co2peaks_yang,interpTime);
interpCo2PeaksCuimei = interp1(timepeaks_cuimei,co2peaks_cuimei,interpTime);

% for each index, select the maximum co2 value
maxCo2Curve = [(interpCo2PeaksYang)',(interpCo2PeaksCuimei)'];
maxCo2Curve = (max(maxCo2Curve'))'; % make maxCurve a signle column

% rename maxCo2Curve as co2envlp
co2envlp    = maxCo2Curve;

% rename interpTime to timepeaks and make it a column vector
timepeaks = interpTime';


% plot
h = figure('pos',[300 600 1500 300]); hold on;
plot(timeline,co2value,'-','LineWidth',1,'color',[0.5 0.5 0.5]);
plot(timeline,criterion * 10 + 1,'-','LineWidth',1,'color',[0.9 0.7 0.0]);
% plot(timepeaks,co2peaks,'bo','LineWidth',1);
plot(timepeaks,co2envlp,'r-','LineWidth',1); 
title('EtCO_2 recording and envelop'); xlabel('Time (sec)'); ylabel('EtCO_2 (mmHg)');
axis([min(timeline) max(timeline) 0 max(co2value)*1.1]);
tightfig;

% save figure and envelop in text
saveas(h,[physiodir filesep txtname_co2_raw(1:end-4)],'fig');
dlmwrite([physiodir filesep txtname_co2_raw(1:end-4) '_envelop_total.txt'], [timepeaks co2envlp]);


%% wrap up
txtname_etco2   = [physiodir filesep txtname_co2_raw(1:end-4) '_envelop_total.txt'];


