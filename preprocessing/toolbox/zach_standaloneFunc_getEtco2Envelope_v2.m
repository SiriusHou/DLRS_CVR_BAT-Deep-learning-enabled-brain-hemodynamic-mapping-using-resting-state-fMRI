function [txtname_etco2, timepeaks, co2envlp] = ...
    zach_standaloneFunc_getEtco2Envelope_v2(physiodir,...
    txtname_co2_raw, sample_rate_co2, envelope_interp_rate, select,...
    varargin)
% Reads co2 trace data and finds the End-Tidal CO2 envelope. Writes this
% data into a text file in physiodir with "_envelop_total" postfixed to the
% original file name (txtname_co2_raw)

% INPUT OVERVIEW:
% physiodir             - directory of the co2 trace
% txtname_co2_raw       - name of the co2 trace text file
% sample_rate_co2       - sampling frequency (in Hz) of the co2 trace data
% envelope_interp_rate  - interpolation rate in s of envelope (set as 1)
% select                - control parameter in picking peaks (set as 0.1)
% >>> Optional inputs:
% varargin(1)           - start time (in s) of the co2 trace (omit all data before)
% varargin(2)           - end time (in s) of the co2 trace (omit all data after)

% (e.g.)
% physiodir = 'C:\Users\Zach\workspace';
% txtname_co2_raw = 'co2_recording.txt';
% sample_rate_co2 = 100;
% envelope_interp_rate = 1;
% select = 0.1;
% varargin(1) = 120;
% varargin(2) = 250;

% OUTPUT OVERVIEW:
% txtname_etco2         - path to the text file containing the co2 envelope
% timepeaks             - timestamps of the co2 envelope (s)
% co2envlp              - cooresponding values of the co2 envelope (mmHg)

disp('Finding EtCO2 Envelope...')
%% inputs
fn_co2              = [physiodir filesep txtname_co2_raw];


%% deal with etco2 recording
% load text file
fid         = fopen(fn_co2);
data_co2    = textscan(fid,'%f');
fclose(fid);
vec_co2     = data_co2{1};
timep_etco2 = ( 0:length(vec_co2)-1 )' / sample_rate_co2;

% chop data if start/end times are given
startTime = 0;
endTime = timep_etco2(end);
if length(varargin) > 0
    startTime = varargin{1};
end
if length(varargin) > 1
    endTime = varargin{2};
end
killIndex = find((timep_etco2<=startTime)|(timep_etco2>=endTime));
timep_etco2(killIndex) = [];
vec_co2(killIndex) = [];

% initial smoothing
value_etco2 = smooth(vec_co2,21);
curve_etco2 = [timep_etco2,value_etco2];


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

% %%
% figure;
% plot(timeline,co2value,'-','LineWidth',1,'color',[0.5 0.5 0.5]);
% hold on
% plot(timepeaks_cuimei,co2peaks_cuimei,'b*')
% plot(timepeaks_yang,co2peaks_yang,'r*')
% legend('trace','cuimei','yang')
% %%

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

% %%
% figure;
% plot(timeline,co2value,'-','LineWidth',1,'color',[0.5 0.5 0.5]);
% hold on
% plot(interpTime,interpCo2PeaksCuimei,'bo')
% plot(interpTime,interpCo2PeaksYang,'ro')
% legend('trace','cuimei','yang')
% %%


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




%% accompanying functions

function output = zach_cuimeiEtco2Alg(timeData,co2Data)
% takes a column of timestamps and a column of co2 levels and finds the
% end-tidal co2 peaks. This is done by smoothing the co2 trace and then
% using the crosses between the smoothed trace and the actual trace to
% search for local maxima.

    % get CO2t and COt out and get rid of the0
    CO2t=timeData;
    CO2=co2Data;
%     for idx=2:size(CO2t,1)
%         if (CO2t(idx)<0.001)
%             CO2t=CO2t(1:idx-1);
%             CO2=CO2(1:idx-1);
%             break;
%         end
%     end
    SmoothCO2=smooth(CO2,100); %10 is the smooth span need to otimize
    
    %looking for cross
    imax=0;
    for idx=2:size(CO2,1)-1
        if ((SmoothCO2(idx)-CO2(idx))>0)&((SmoothCO2(idx-1)-CO2(idx-1))<0)
            imax=imax+1;
            areaidx(imax)=idx;
            if (imax==1)

                 [peakCO2(imax,2), maxidx]=max(CO2(1:idx));
                 tmpCOt=CO2t(1:idx);
                  peakCO2(imax,1)=tmpCOt(maxidx);
            else
                 [peakCO2(imax,2), maxidx] =max(CO2(areaidx(imax-1):idx));
                    tmpCOt=CO2t(areaidx(imax-1):idx);
                  peakCO2(imax,1)=tmpCOt(maxidx);
            end
        end
    end
    
    output = peakCO2;

end

function hfig = tightfig(hfig)
% tightfig: Alters a figure so that it has the minimum size necessary to
% enclose all axes in the figure without excess space around them.
% 
% Note that tightfig will expand the figure to completely encompass all
% axes if necessary. If any 3D axes are present which have been zoomed,
% tightfig will produce an error, as these cannot easily be dealt with.
% 
% Input
%
% hfig - handle to figure, if not supplied, the current figure will be used
%   instead.
%
%

    if nargin == 0
        hfig = gcf;
    end

    % There can be an issue with tightfig when the user has been modifying
    % the contnts manually, the code below is an attempt to resolve this,
    % but it has not yet been satisfactorily fixed
%     origwindowstyle = get(hfig, 'WindowStyle');
    set(hfig, 'WindowStyle', 'normal');
    
    % 1 point is 0.3528 mm for future use

    % get all the axes handles note this will also fetch legends and
    % colorbars as well
    hax = findall(hfig, 'type', 'axes');
    % TODO: fix for modern matlab, colorbars and legends are no longer axes
    hcbar = findall(hfig, 'type', 'colorbar');
    hleg = findall(hfig, 'type', 'legend');
    
    % get the original axes units, so we can change and reset these again
    % later
    origaxunits = get(hax, 'Units');
    
    % change the axes units to cm
    set(hax, 'Units', 'centimeters');
    
    pos = [];
    ti = [];
    
    % get various position parameters of the axes
    if numel(hax) > 1
%         fsize = cell2mat(get(hax, 'FontSize'));
        ti = cell2mat(get(hax,'TightInset'));
        pos = [pos; cell2mat(get(hax, 'Position')) ];
    else
%         fsize = get(hax, 'FontSize');
        ti = get(hax,'TightInset');
        pos = [pos; get(hax, 'Position') ];
    end
    
    if ~isempty (hcbar)
        
        set(hcbar, 'Units', 'centimeters');
        
        % colorbars do not have tightinset property
        for cbind = 1:numel(hcbar)
            %         fsize = cell2mat(get(hax, 'FontSize'));
            [cbarpos, cbarti] = colorbarpos (hcbar);

            pos = [pos; cbarpos];
            ti = [ti; cbarti];
        end
    end
    
    if ~isempty (hleg)
        
        set(hleg, 'Units', 'centimeters');
        
        % legends do not have tightinset property
        if numel(hleg) > 1
            %         fsize = cell2mat(get(hax, 'FontSize'));
            pos = [pos; cell2mat(get(hleg, 'Position')) ];
        else
            %         fsize = get(hax, 'FontSize');
            pos = [pos; get(hleg, 'Position') ];
        end
        ti = [ti; repmat([0,0,0,0], numel(hleg), 1); ];
    end
    
    % ensure very tiny border so outer box always appears
    ti(ti < 0.1) = 0.15;
    
    % we will check if any 3d axes are zoomed, to do this we will check if
    % they are not being viewed in any of the 2d directions
    views2d = [0,90; 0,0; 90,0];
    
    for i = 1:numel(hax)
        
        set(hax(i), 'LooseInset', ti(i,:));
%         set(hax(i), 'LooseInset', [0,0,0,0]);
        
        % get the current viewing angle of the axes
        [az,el] = view(hax(i));
        
        % determine if the axes are zoomed
        iszoomed = strcmp(get(hax(i), 'CameraViewAngleMode'), 'manual');
        
        % test if we are viewing in 2d mode or a 3d view
        is2d = all(bsxfun(@eq, [az,el], views2d), 2);
               
        if iszoomed && ~any(is2d)
           error('TIGHTFIG:haszoomed3d', 'Cannot make figures containing zoomed 3D axes tight.') 
        end
        
    end
    
    % we will move all the axes down and to the left by the amount
    % necessary to just show the bottom and leftmost axes and labels etc.
    moveleft = min(pos(:,1) - ti(:,1));
    
    movedown = min(pos(:,2) - ti(:,2));
    
    % we will also alter the height and width of the figure to just
    % encompass the topmost and rightmost axes and lables
    figwidth = max(pos(:,1) + pos(:,3) + ti(:,3) - moveleft);
    
    figheight = max(pos(:,2) + pos(:,4) + ti(:,4) - movedown);
    
    % move all the axes
    for i = 1:numel(hax)
        
        set(hax(i), 'Position', [pos(i,1:2) - [moveleft,movedown], pos(i,3:4)]);
        
    end
    
    for i = 1:numel(hcbar)
        
        set(hcbar(i), 'Position', [pos(i+numel(hax),1:2) - [moveleft,movedown], pos(i+numel(hax),3:4)]);
        
    end
    
    for i = 1:numel(hleg)
        
        set(hleg(i), 'Position', [pos(i+numel(hax)+numel(hcbar),1:2) - [moveleft,movedown], pos(i+numel(hax)+numel(hcbar),3:4)]);
        
    end
    
    origfigunits = get(hfig, 'Units');
    
    set(hfig, 'Units', 'centimeters');
    
    % change the size of the figure
    figpos = get(hfig, 'Position');
    
    set(hfig, 'Position', [figpos(1), figpos(2), figwidth, figheight]);
    
    % change the size of the paper
    set(hfig, 'PaperUnits','centimeters');
    set(hfig, 'PaperSize', [figwidth, figheight]);
    set(hfig, 'PaperPositionMode', 'manual');
    set(hfig, 'PaperPosition',[0 0 figwidth figheight]);    
    
    % reset to original units for axes and figure 
    if ~iscell(origaxunits)
        origaxunits = {origaxunits};
    end

    for i = 1:numel(hax)
        set(hax(i), 'Units', origaxunits{i});
    end

    set(hfig, 'Units', origfigunits);
    
%      set(hfig, 'WindowStyle', origwindowstyle);
     
end

function [pos, ti] = colorbarpos (hcbar)

    % 1 point is 0.3528 mm
    
    pos = hcbar.Position;
    ti = [0,0,0,0];
    
    if ~isempty (strfind (hcbar.Location, 'outside'))

        if strcmp (hcbar.AxisLocation, 'out')
            
            tlabels = hcbar.TickLabels;
            
            fsize = hcbar.FontSize;
            
            switch hcbar.Location
                
                case 'northoutside'
                    
                    % make exta space a little more than the font size/height
                    ticklablespace_cm = 1.1 * (0.3528/10) * fsize;
                    
                    ti(4) = ti(4) + ticklablespace_cm;
                    
                case 'eastoutside'
                    
                    maxlabellen = max ( cellfun (@numel, tlabels, 'UniformOutput', true) );
            
                    % 0.62 factor is arbitrary and added because we don't
                    % know the width of every character in the label, the
                    % fsize refers to the height of the font
                    ticklablespace_cm = (0.3528/10) * fsize * maxlabellen * 0.62;

                    ti(3) = ti(3) + ticklablespace_cm;
                    
                case 'southoutside'
                    
                    % make exta space a little more than the font size/height
                    ticklablespace_cm = 1.1 * (0.3528/10) * fsize;

                    ti(2) = ti(2) + ticklablespace_cm;
                    
                case 'westoutside'
                    
                    maxlabellen = max ( cellfun (@numel, tlabels, 'UniformOutput', true) );
            
                    % 0.62 factor is arbitrary and added because we don't
                    % know the width of every character in the label, the
                    % fsize refers to the height of the font
                    ticklablespace_cm = (0.3528/10) * fsize * maxlabellen * 0.62;

                    ti(1) = ti(1) + ticklablespace_cm;
                    
            end
            
        end
        
    end

end

function [peaks,criterion] = pickpeaks(V,select,display)
% -------------------------------------------------------------
% Scale-space peak picking
% ------------------------
% This function looks for peaks in the data using scale-space theory. 
% 
% input : 
%   * V : data, a vector
%   * select : either: 
%       - select >1 : the number of peaks to detect
%       - 0<select<1 : the threshold to apply for finding peaks 
%         the closer to 1, the less peaks, the closer to 0, the more peaks
%   * display : whether or not to display a figure for the results. 0 by
%               default
%   * ... and that's all ! that's the cool thing about the algorithm =)
% 
% outputs :
%   * peaks : indices of the peaks
%   * criterion : the value of the computed criterion. Same
%                 length as V and giving for each point a high value if
%                 this point is likely to be a peak
% 
% The algorithm goes as follows:
% 1°) set a smoothing horizon, initially 1;
% 2°) smooth the data using this horizon
% 3°) find local extrema of this smoothed data
% 4°) for each of these local extrema, link it to a local extremum found in
%     the last iteration. (initially just keep them all) and increment the 
%     corresponding criterion using current scale. The
%     rationale is that a trajectory surviving such smoothing is an important
%     peak
% 5°) Iterate to step 2°) using a larger horizon.
% 
% At the end, we keep the points with the largest criterion as peaks.
% I don't know if that kind of algorithm has already been published
% somewhere, I coded it myself and it works pretty nice, so.. enjoy !
% If you find it useful, please mention it in your studies by referencing
% the following research report:
%
%@techreport{liutkus:hal-01103123,
%  TITLE = {{Scale-Space Peak Picking}},
%  AUTHOR = {Liutkus, Antoine},
%  URL = {https://hal.inria.fr/hal-01103123},
%  TYPE = {Research Report},
%  INSTITUTION = {{Inria Nancy - Grand Est (Villers-l{\`e}s-Nancy, France)}},
%  YEAR = {2015},
%  MONTH = Jan,
%  KEYWORDS = { scale-space ; peak detection},
%  HAL_ID = {hal-01103123},
%  HAL_VERSION = {v1},
%}
%
% 
% running time should be decent, although intrinsically higher than 
% findpeaks. For vectors of length up to, say, 10 000, it should be nice. 
% Above, it may be worth it though.
% ---------------------------------------------------------------------
% Copyright (C) 2015, Inria, Antoine Liutkus
%
%Redistribution and use in source and binary forms, with or without
%modification, are permitted provided that the following conditions are met:
%    * Redistributions of source code must retain the above copyright
%      notice, this list of conditions and the following disclaimer.
%    * Redistributions in binary form must reproduce the above copyright
%      notice, this list of conditions and the following disclaimer in the
%       documentation and/or other materials provided with the distribution.
%     * Neither the name of Inria nor the names of its contributors may 
%       be used to endorse or promote products derived from this software 
%       without specific prior written permission.
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
% ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
% WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL INRIA BE LIABLE FOR ANY
% DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
% (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
% LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
% ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
% (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

%data is a vector
V = V(:)-min((V(:)));

%input parsin
if nargin < 3
    display=0;
end
if nargin < 2
    select= 0;
end

n = length(V);

%definition of local variables
buffer = zeros(n,1);
criterion = zeros(n,1);
if select < 1
    minDist = n/20;
else
    minDist = n/select;
end
%horizons = round(linspace(1,ceil(n/20),50));
horizons = unique(round(logspace(0,2,50)/100*ceil(n/20)));

%horizons=1:2:50;
Vorig = V;

% all this tempMat stuff is to avoid calling findpeaks which is horribly
% slow for our purpose
tempMat = zeros(n,3);
tempMat(1,1)=inf;
tempMat(end,3)=inf;

% loop over scales
for is=1:length(horizons)
    
    %sooth data, using fft-based convolution with a half sinusoid
    horizon = horizons(is);
    if horizon > 1
        w=max(eps,sin(2*pi*(0:(horizon-1))/2/(horizon-1)));
        w=w/sum(w);    
        %V=conv(V(:),w(:),'same');
        V = real(ifft(fft(V(:),n+horizon).*fft(w(:),n+horizon)));
        V = V(1+floor(horizon/2):end-ceil(horizon/2));
    end

    %find local maxima
    tempMat(2:end,1) = V(1:end-1);
    tempMat(:,2) = V(:);
    tempMat(1:end-1,3) = V(2:end);
    [useless,posMax] =max(tempMat,[],2);
    I = find(posMax==2);
    I = I(:)';
    
    %initialize buffer
    newBuffer = zeros(size(buffer));
    
    if is == 1
        % if first iteration, keep all local maxima
        newBuffer(I) = Vorig(I);
    else    
        old = find(buffer);
        old = old(:)';
        if isempty(old)
            continue;
        end
        
        %Now, for each element of I, find the closest element in
        %old along with its distance. The few nice lines below were
        %written by Roger Stafford in a forum post available here:
        %http://www.mathworks.fr/matlabcentral/newsreader/view_thread/24387
        [c,p] = sort(old);
        [useless,ic] = histc(I,[-inf,(c(1:end-1)+c(2:end))/2,inf]);
        iOld = p(ic);
        d = abs(I-old(iOld));
 
        %done, now select only those that are sufficiently close
        neighbours = iOld(d<minDist);

        if ~isempty(neighbours)
            newBuffer(old(neighbours)) = V(old(neighbours))*is^2;
        end
    end
    %update stuff
    buffer = newBuffer;
    criterion = criterion + newBuffer;
end

%normalize criterion
criterion = criterion/max(criterion);

%find peaks based on criterion
if select<1
    peaks = find(criterion>select);
else
%     sorted = find(criterion>1E-3);
%     [~,order] = sort(criterion(sorted),'descend');
%     peaks = sorted(order(1:min(length(sorted),select)));
    [useless,order] = sort(criterion,'descend');
    peaks = order(1:select);
end

if display
    %display
    clf
    plot(Vorig,'LineWidth',2);
    hold on
    plot(criterion*max(Vorig),'r');
    hold on
    plot(peaks,Vorig(peaks),'ro','MarkerSize',10,'LineWidth',2)
    grid on
    title('Scale-space peak detection','FontSize',16);
    legend('data','computed criterion','selected peaks');
end
end

end
