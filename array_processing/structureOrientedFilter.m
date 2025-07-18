function DataStruct = structureOrientedFilter(DataStruct, param)
% RADONTRANSFORM  Perform Radon Transform-based array processing on DataStruct.
%
% Usage:
%   DataStruct = radonTransform(DataStruct, param)
%
% Inputs:
%   DataStruct : struct array with fields:
%       .EventInfo.evid         - Event ID (string)
%       .Waveforms.dataProcessed - Processed waveform data [Nt x 3]
%       .RF.ittime              - (not always used here, but often stored in RF)
%       .TravelInfo.distDeg     - Epicentral distance in degrees
%       .ProcHistory            - Cell array to record logs (optional)
%
%   param : struct with (optional) fields:
%       .lows      (default: 0.1)    - Low corner freq for bandpass (Hz)
%       .highs     (default: 1.2)    - High corner freq for bandpass (Hz)
%       .pmax      (default: 0.05)   - Max slowness (s/km)
%       .pmin      (default: -0.05)  - Min slowness (s/km)
%       .minTraces (default: 60)     - Min number of traces per event
%       .N1        (default: 30)     - Iterations for CG solver
%       .N2        (default: 1)      - Sparse solution weight
%       .plotRadon (default: false)  - Flag to plot Radon results
%
% Output:
%   DataStruct : The updated struct array. For each trace that belongs
%                to an event with enough traces, the field:
%       DataStruct(n).Waveforms.dataStructureOrientedFiltered
%                is populated with [Nt x 3] radial/transverse (or vertical)
%                waveforms after Radon-based noise removal or separation.
%
% Author:  Yunfeng Chen (Refined by ChatGPT)
% Date   : Jan. 12, 2025

%% 1. Parameter Handling & Validation
if nargin < 2 || isempty(param)
    param = struct();
end

% Define default parameters if missing
if ~isfield(param, 'lows'),       param.lows       = 0.1;   end
if ~isfield(param, 'highs'),      param.highs      = 1.2;   end
if ~isfield(param, 'pmax'),       param.pmax       = 0.05;  end
if ~isfield(param, 'pmin'),       param.pmin       = -0.05; end
if ~isfield(param, 'minTraces'),  param.minTraces  = 60;    end
if ~isfield(param, 'plotRadon'),  param.plotRadon  = false; end

% Basic field checks (on the first element, assuming consistent struct array)
requiredTopLevel = {'EventInfo','Waveforms','TravelInfo','RF'};
for f = requiredTopLevel
    if ~isfield(DataStruct, f{1})
        error('radonTransform:MissingField',...
            'DataStruct must contain the field %s in each element.', f{1});
    end
end
if ~isfield(DataStruct(1).EventInfo, 'evid')
    error('radonTransform:MissingField',...
        'DataStruct.EventInfo must contain field: evid.');
end
if ~isfield(DataStruct(1).Waveforms, 'dataProcessed')
    error('radonTransform:MissingField',...
        'DataStruct.Waveforms must contain field: dataProcessed.');
end
if ~isfield(DataStruct(1).TravelInfo, 'distDeg')
    error('radonTransform:MissingField',...
        'DataStruct.TravelInfo must contain field: distDeg.');
end
%% 3. Identify Unique Events
eventListStruct = getEvents(DataStruct);  % retrieve unique events from the data
eventIDs = {eventListStruct.evid};

%% 4. Process Each Event
for iEvt = 1:length(eventIDs)
    eventID = eventIDs{iEvt};
    [commonEventGather, matchIndex] = getCommonEventGather(DataStruct, eventID);

    % Check if enough traces are available for Radon transform
    if length(commonEventGather) < param.minTraces
        skipMsg = sprintf('[RadonTransform] Event %s skipped: only %d traces (< %d).',...
            eventID, length(commonEventGather), param.minTraces);
        commonEventGather = appendHistory(commonEventGather, skipMsg);
        DataStruct(matchIndex) = commonEventGather;
        continue;
    end

    % Display progress
    fprintf('Processing event #%d: %s with %d traces.\n', ...
        iEvt, eventID, length(commonEventGather));

    dist = zeros(length(commonEventGather),1);
    for ii = 1:length(commonEventGather)
        deg  = commonEventGather(ii).TravelInfo.distDeg;
        dist(ii) = deg * (2*pi*6371/360);  % degrees -> km
    end

    % Convert offsets to h = distance - min(distance)
    h = dist - min(dist);
    [h,idx] = sort(h);
    commonEventGather = commonEventGather(idx);

    %% 4.1 Extract and Preprocess Waveforms
    [d_z, d_r, d_t, t, dt] = extractAndPreprocessWaveforms(commonEventGather,...
        param.lows, param.highs);

    %% 4.3 Perform Radon Transform on Z
    try
        %Calculate structural dips
        dip = str_dip2d(d_z);

        % Apply structure-oriented smoothing
        %         dp_z = str_pwsmooth_lop2d(d_z, dip, 3, 2, 0.01);
        [n1,n2] = size(d_z);
        type_mf=0;
        ifsmooth=0;
        ndn=n1*n2;
        nds=n1*n2;
        ns=3;
        order=2;
        eps=0.01;
        [~,dp_z]=amf_pwsmooth_lop(dip,[],n1,n2,ns,order,eps,ndn,nds,type_mf,ifsmooth,d_z,[]);
        dp_z = reshape(dp_z,n1,n2);
    catch ME
        warnMsg = sprintf('[RadonTransform] Event %s, Z-component failed: %s', ...
            eventID, ME.message);
        warning(warnMsg);
        commonEventGather = appendHistory(commonEventGather, warnMsg);
        DataStruct(matchIndex) = commonEventGather;
        continue;
    end

    %% 4.4 Perform Radon Transform on R
    try
        dip = str_dip2d(d_r);

        % Apply structure-oriented smoothing
        %         dp_r = str_pwsmooth_lop2d(d_r, dip, 3, 2, 0.01);
        [n1,n2] = size(d_r);
        type_mf=0;
        ifsmooth=0;
        ndn=n1*n2;
        nds=n1*n2;
        ns=3;
        order=2;
        eps=0.01;
        [~,dp_r]=amf_pwsmooth_lop(dip,[],n1,n2,ns,order,eps,ndn,nds,type_mf,ifsmooth,d_r,[]);
        dp_r = reshape(dp_r,n1,n2);
    catch ME
        warnMsg = sprintf('[RadonTransform] Event %s, R-component failed: %s', ...
            eventID, ME.message);
        warning(warnMsg);
        commonEventGather = appendHistory(commonEventGather, warnMsg);
        DataStruct(matchIndex) = commonEventGather;
        continue;
    end

    %% 4.6 Store Radon-Filtered Waveforms
    % If dp_t is not computed, you can fill T with zeros or skip entirely
    for n = 1:length(commonEventGather)
        % Initialize dataStructureOrientedFiltered if needed
        ntData = size(d_z,1);  % or length(t)
        commonEventGather(n).Waveforms.dataStructureOrientedFiltered = zeros(ntData, 3);

        % T = dp_t(:,n), or zeros if T not used
        commonEventGather(n).Waveforms.dataStructureOrientedFiltered(:,1) = 0;
        % R
        commonEventGather(n).Waveforms.dataStructureOrientedFiltered(:,2) = dp_r(:,n);
        % Z
        commonEventGather(n).Waveforms.dataStructureOrientedFiltered(:,3) = dp_z(:,n);
    end

    % Log the successful Radon transform
    successMsg = sprintf('[RadonTransform] Event %s: Radon transform applied (%d traces).',...
        eventID, length(commonEventGather));
    commonEventGather = appendHistory(commonEventGather, successMsg);

    % Place updated event gather back into DataStruct
    DataStruct(matchIndex) = commonEventGather;

    %% 4.7 Optional: Plot Radon Results
    if param.plotRadon
        plotRadonResults(d_z, d_r, dp_z, dp_r, h, t, eventID);
    end
end

end % end of radonTransform main function

%% ------------------------- LOCAL SUBFUNCTIONS -------------------------%%

function [d_z, d_r, d_t, t, dt] = extractAndPreprocessWaveforms(commonEventGather, lows, highs)
% EXTRACTANDPREPROCESSWAVEFORMS  Retrieve and bandpass-filter Z/R/T from
%                                each trace in the gather.
%
% Inputs:
%   commonEventGather : array of DataStruct elements for the same event
%   lows, highs       : bandpass corners
%
% Outputs:
%   d_z, d_r, d_t : [Nt x Ntraces] arrays for each component
%   t            : time vector [Nt x 1] from the first trace
%   dt           : sample spacing (sec)

nTraces = length(commonEventGather);
dataSample = commonEventGather(1).Waveforms.dataProcessed;
[nt, ~] = size(dataSample);

% Prepare storage
d_z = zeros(nt, nTraces);
d_r = zeros(nt, nTraces);
d_t = zeros(nt, nTraces);

% Time step from the first trace (assumes all are the same)
tVec = commonEventGather(1).TimeAxis.t_resample;
dt = tVec(2) - tVec(1);

% Process each trace
for iTr = 1:nTraces
    dataProc = commonEventGather(iTr).Waveforms.dataProcessed;  % [Nt x 3]
    % Components: T=1, R=2, Z=3
    tmpZ = dataProc(:,3);
    tmpR = dataProc(:,2);
    tmpT = dataProc(:,1);

    % Taper (5 sec at start/end) - adapt function arguments as needed
    tmpZ = taper(tmpZ, 5, 5);
    tmpR = taper(tmpR, 5, 5);
    tmpT = taper(tmpT, 5, 5);

    % Bandpass filter
    tmpZ = bandpassSeis(tmpZ, dt, lows, highs, 3);
    tmpR = bandpassSeis(tmpR, dt, lows, highs, 3);
    tmpT = bandpassSeis(tmpT, dt, lows, highs, 3);

    % Store in output arrays
    d_z(:, iTr) = tmpZ;
    d_r(:, iTr) = tmpR;
    d_t(:, iTr) = tmpT;
end

t = tVec;  % consistent for all traces
end

function gatherStruct = appendHistory(gatherStruct, msg)
% APPENDHISTORY  Append a message to each element's ProcHistory in gatherStruct.
% If ProcHistory is empty or missing, initialize it.

for ii = 1:length(gatherStruct)
    if ~isfield(gatherStruct(ii), 'ProcHistory') || isempty(gatherStruct(ii).ProcHistory)
        gatherStruct(ii).ProcHistory = {msg};
    else
        gatherStruct(ii).ProcHistory{end+1} = msg;
    end
end
end

function plotRadonResults(d_z, d_r, dp_z, dp_r, h, t, eventID)
% PLOTRADONRESULTS  Visualize raw vs. Radon-separated data for Z and R.
%
% Inputs:
%   d_z, dp_z : raw vs. Radon-transformed data for Z ( [Nt x Ntraces] )
%   d_r, dp_r : raw vs. Radon-transformed data for R
%   h         : distance offsets (km) for each trace
%   t         : time vector (s)
%   eventID   : string identifier for the event

% Estimate amplitude for color scaling
cmax_z = 3 * rms(dp_z(:));
cmax_r = 3 * rms(dp_r(:));

figure('Name', sprintf('Structure Oriented Filter Results for Event %s', eventID), ...
    'Position', [100, 100, 1200, 800], 'Color', 'w');

% Raw Z
subplot(2,3,1);
imagesc(h, t, d_z);
colormap(seismic(1));
caxis([-cmax_z cmax_z]);
xlabel('Distance (km)'); ylabel('Time (s)');
title('Raw Z'); set(gca, 'FontSize', 12);

% Radon Z
subplot(2,3,2);
imagesc(h, t, dp_z);
colormap(seismic(1));
caxis([-cmax_z cmax_z]);
xlabel('Distance (km)'); ylabel('Time (s)');
title('Radon Z'); set(gca, 'FontSize', 12);

% Z difference
subplot(2,3,3);
imagesc(h, t, d_z - dp_z);
colormap(seismic(1));
caxis([-cmax_z cmax_z]);
xlabel('Distance (km)'); ylabel('Time (s)');
title('Removed Noise (Z)'); set(gca, 'FontSize', 12);

% Raw R
subplot(2,3,4);
imagesc(h, t, d_r);
colormap(seismic(1));
caxis([-cmax_r cmax_r]);
xlabel('Distance (km)'); ylabel('Time (s)');
title('Raw R'); set(gca, 'FontSize', 12);

% Radon R
subplot(2,3,5);
imagesc(h, t, dp_r);
colormap(seismic(1));
caxis([-cmax_r cmax_r]);
xlabel('Distance (km)'); ylabel('Time (s)');
title('Radon R'); set(gca, 'FontSize', 12);

% R difference
subplot(2,3,6);
imagesc(h, t, d_r - dp_r);
colormap(seismic(1));
caxis([-cmax_r cmax_r]);
xlabel('Distance (km)'); ylabel('Time (s)');
title('Removed Noise (R)'); set(gca, 'FontSize', 12);
end
