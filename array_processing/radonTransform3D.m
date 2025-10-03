function DataStruct = radonTransform3D(DataStruct, gridStruct, param)
% RADONTRANSFORM3D  Perform 3D Radon Transform-based array processing on DataStruct.
%
% Usage:
%   DataStruct = radonTransform3D(DataStruct, gridStruct, param)
%
% Inputs:
%   DataStruct : struct array with fields:
%       .EventInfo.evid         - Event ID (string)
%       .Waveforms.dataProcessed - Processed waveform data [Nt x 3]
%       .RF.ittime              - (not always used here, but often stored in RF)
%       .TravelInfo.distDeg     - Epicentral distance in degrees
%       .ProcHistory            - Cell array to record logs (optional)
%   gridStruct : struct with grid information containing station coordinates
%   param : struct with (optional) fields:
%       .lows      (default: 0.1)    - Low corner freq for bandpass (Hz)
%       .highs     (default: 1.2)    - High corner freq for bandpass (Hz)
%       .pxmax     (default: 0.05)   - Max x-slowness (s/km)
%       .pxmin     (default: -0.05)  - Min x-slowness (s/km)
%       .pymax     (default: 0.05)   - Max y-slowness (s/km)
%       .pymin     (default: -0.05)  - Min y-slowness (s/km)
%       .minTraces (default: 60)     - Min number of traces per event
%       .N1        (default: 30)     - Iterations for CG solver
%       .N2        (default: 1)      - Sparse solution weight
%       .plotRadon (default: false)  - Flag to plot Radon results
%       .order     (default: postdecon) - Apply of filtering before or after deconvolution
%       .type      (default: 1)      - Radon type (1=linear, 2=parabolic, 3=hyperbolic)
%
% Output:
%   DataStruct : The updated struct array. For each trace that belongs
%                to an event with enough traces, the field:
%       DataStruct(n).Waveforms.dataRadonFiltered
%                is populated with [Nt x 3] radial/transverse (or vertical)
%                waveforms after Radon-based noise removal or separation.
%
% Author:  Based on 2D version by Yunfeng Chen, extended to 3D
% Date   : Oct. 3, 2025

%% 1. Parameter Handling & Validation
if nargin < 3 || isempty(param)
    param = struct();
end

% Define default parameters if missing
if ~isfield(param, 'lows'),       param.lows       = 0.1;   end
if ~isfield(param, 'highs'),      param.highs      = 1.2;   end
if ~isfield(param, 'pxmax'),      param.pxmax      = 0.05;  end
if ~isfield(param, 'pxmin'),      param.pxmin      = -0.05; end
if ~isfield(param, 'pymax'),      param.pymax      = 0.05;  end
if ~isfield(param, 'pymin'),      param.pymin      = -0.05; end
if ~isfield(param, 'minTraces'),  param.minTraces  = 60;    end
if ~isfield(param, 'N1'),         param.N1         = 30;    end
if ~isfield(param, 'N2'),         param.N2         = 1;     end
if ~isfield(param, 'plotRadon'),  param.plotRadon  = false; end
if ~isfield(param, 'order'),      param.order      = 'postdecon'; end
if ~isfield(param, 'type'),       param.type       = 1;     end

% Basic field checks (on the first element, assuming consistent struct array)
requiredTopLevel = {'EventInfo','Waveforms','TravelInfo','RF'};
for f = requiredTopLevel
    if ~isfield(DataStruct, f{1})
        error('radonTransform3D:MissingField',...
            'DataStruct must contain the field %s in each element.', f{1});
    end
end
if ~isfield(DataStruct(1).EventInfo, 'evid')
    error('radonTransform3D:MissingField',...
        'DataStruct.EventInfo must contain field: evid.');
end
if ~isfield(DataStruct(1).Waveforms, 'dataProcessed')
    error('radonTransform3D:MissingField',...
        'DataStruct.Waveforms must contain field: dataProcessed.');
end
if ~isfield(DataStruct(1).TravelInfo, 'distDeg')
    error('radonTransform3D:MissingField',...
        'DataStruct.TravelInfo must contain field: distDeg.');
end

%% 2. 3D Radon Parameter Setup
% Define slowness axes for x and y directions (s/km)
px = linspace(param.pxmin, param.pxmax, 25);
py = linspace(param.pymin, param.pymax, 25);
npx = length(px);
npy = length(py);

% Create Param structure to pass to yc_pcg() or radon3d_op()
Param.px = px;
Param.py = py;
Param.N1 = param.N1;
Param.N2 = param.N2;
Param.type = param.type;

%% 3. Identify Unique Events
eventListStruct = getEvents(DataStruct);  % retrieve unique events from the data
eventIDs = {eventListStruct.evid};

%% 4. Process Each Event
for iEvt = 1:length(eventIDs)
    eventID = eventIDs{iEvt};
    [commonEventGather, matchIndex] = getCommonEventGather(DataStruct, eventID);

    % Check if enough traces are available for Radon transform
    if length(commonEventGather) < param.minTraces
        skipMsg = sprintf('[RadonTransform3D] Event %s skipped: only %d traces (< %d).',...
            eventID, length(commonEventGather), param.minTraces);
        commonEventGather = appendHistory3D(commonEventGather, skipMsg);
        DataStruct(matchIndex) = commonEventGather;
        continue;
    end

    % Display progress
    fprintf('Processing event #%d: %s with %d traces.\n', ...
        iEvt, eventID, length(commonEventGather));

    %% 4.1 Get 2D Station Coordinates
    stationList = getStations(DataStruct);
    stlo = [stationList.stlo]';  % station longitude
    stla = [stationList.stla]';  % station latitude
    
    % Convert to projected coordinates (2D)
    [hx, hy] = latlonToProjectedCoords(stlo, stla, gridStruct);
    
    % Sort by x coordinate for consistent ordering
    [hx_sorted, idx] = sort(hx);
    hy_sorted = hy(idx);
    commonEventGather = commonEventGather(idx);

    %% 4.2 Extract and Preprocess Waveforms
    [d_z, d_r, d_t, t, dt] = extractAndPreprocessWaveforms3D(commonEventGather,...
        param.lows, param.highs);
    trace_energy = rms(d_z);
    % remove traces with anomalous amplitude
    remove_idx = trace_energy > mean(trace_energy) + 2*std(trace_energy);
    d_z(:,remove_idx) = 0;
    d_r(:,remove_idx) = 0;

    switch param.order
        case 'predecon'
            %% Set 3D radon parameters
            Param.hx = hx_sorted;     % x-offsets
            Param.hy = hy_sorted;     % y-offsets  
            Param.nt = length(t);
            Param.dt = dt;

            % Preallocate transform model "ma" based on type
            if Param.type == 3  % Hyperbolic
                nv = 25;  % default number of velocities
                Param.v = linspace(1, 8, nv);  % velocity range 1-8 km/s
                ma = zeros(Param.nt, nv);
            else  % Linear or Parabolic
                ma = zeros(Param.nt, npx, npy);
            end

            %% 4.3 Perform 3D Radon Transform on Z
            try
                % Initialize the model with zeros
                mi_z = yc_pcg(@radon3d_op, Param, d_z, ma, param.N1, param.N2, 1);
                dp_z = radon3d_op(mi_z, Param, 1);  % forward modeling from the found model

            catch ME
                warnMsg = sprintf('[RadonTransform3D] Event %s, Z-component failed: %s', ...
                    eventID, ME.message);
                warning(warnMsg);
                commonEventGather = appendHistory3D(commonEventGather, warnMsg);
                DataStruct(matchIndex) = commonEventGather;
                continue;
            end

            %% 4.4 Perform 3D Radon Transform on R
            try
                mi_r = yc_pcg(@radon3d_op, Param, d_r, ma, param.N1, param.N2, 1);
                dp_r = radon3d_op(mi_r, Param, 1);

            catch ME
                warnMsg = sprintf('[RadonTransform3D] Event %s, R-component failed: %s', ...
                    eventID, ME.message);
                warning(warnMsg);
                commonEventGather = appendHistory3D(commonEventGather, warnMsg);
                DataStruct(matchIndex) = commonEventGather;
                continue;
            end

            %% 4.5 Store Radon-Filtered Waveforms
            for n = 1:length(commonEventGather)
                % Initialize dataRadonFiltered if needed
                ntData = size(d_z,1);
                commonEventGather(n).Waveforms.dataRadonFiltered = zeros(ntData, 3);

                % T = zeros (not processed), R and Z from Radon transform
                commonEventGather(n).Waveforms.dataRadonFiltered(:,1) = 0;
                commonEventGather(n).Waveforms.dataRadonFiltered(:,2) = dp_r(:,n);
                commonEventGather(n).Waveforms.dataRadonFiltered(:,3) = dp_z(:,n);
            end

            % Log the successful Radon transform
            successMsg = sprintf('[RadonTransform3D] Event %s: 3D Radon transform applied (%d traces).',...
                eventID, length(commonEventGather));
            commonEventGather = appendHistory3D(commonEventGather, successMsg);

            % Place updated event gather back into DataStruct
            DataStruct(matchIndex) = commonEventGather;

            %% 4.6 Optional: Plot Radon Results
            if param.plotRadon
                plotRadonResults3D(d_z, d_r, dp_z, dp_r, hx_sorted, hy_sorted, t, eventID);
            end

        case 'postdecon'
            %% Perform 3D Radon Transform on RF
            itrCell = {commonEventGather.RF};
            d = cell2mat(cellfun(@(rf) rf.itr, itrCell,'UniformOutput', false));
            d(:,remove_idx) = 0;
            t = commonEventGather(1).RF.ittime;
            
            % Set 3D radon parameters
            Param.hx = hx_sorted;
            Param.hy = hy_sorted;
            Param.nt = length(t);
            Param.dt = dt;

            % Preallocate transform model "ma" based on type
            if Param.type == 3  % Hyperbolic
                nv = 25;  % default number of velocities
                Param.v = linspace(1, 8, nv);  % velocity range 1-8 km/s
                ma = zeros(Param.nt, nv);
            else  % Linear or Parabolic
                ma = zeros(Param.nt, npx, npy);
            end

            try
                % Initialize the model with zeros
                mi = yc_pcg(@radon3d_op, Param, d, ma, param.N1, param.N2, 1);
                dp = radon3d_op(mi, Param, 1);  % forward modeling from the found model

            catch ME
                warnMsg = sprintf('[RadonTransform3D] Event %s, RF failed: %s', ...
                    eventID, ME.message);
                warning(warnMsg);
                commonEventGather = appendHistory3D(commonEventGather, warnMsg);
                DataStruct(matchIndex) = commonEventGather;
                continue;
            end

           for n = 1:length(commonEventGather)
                % save RF
                commonEventGather(n).RF.itr = dp(:,n);
            end

            % Log the successful Radon transform
            successMsg = sprintf('[RadonTransform3D] Event %s: 3D Radon transform applied (%d traces).',...
                eventID, length(commonEventGather));
            commonEventGather = appendHistory3D(commonEventGather, successMsg);

            % Place updated event gather back into DataStruct
            DataStruct(matchIndex) = commonEventGather;
    end
end

end % end of radonTransform3D main function

%% ------------------------- LOCAL SUBFUNCTIONS -------------------------%%

function [d_z, d_r, d_t, t, dt] = extractAndPreprocessWaveforms3D(commonEventGather, lows, highs)
% EXTRACTANDPREPROCESSWAVEFORMS3D  Retrieve and bandpass-filter Z/R/T from
%                                  each trace in the gather for 3D processing.
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

    % Taper (5 sec at start/end)
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

function gatherStruct = appendHistory3D(gatherStruct, msg)
% APPENDHISTORY3D  Append a message to each element's ProcHistory in gatherStruct.
% If ProcHistory is empty or missing, initialize it.

for ii = 1:length(gatherStruct)
    if ~isfield(gatherStruct(ii), 'ProcHistory') || isempty(gatherStruct(ii).ProcHistory)
        gatherStruct(ii).ProcHistory = {msg};
    else
        gatherStruct(ii).ProcHistory{end+1} = msg;
    end
end
end

function plotRadonResults3D(d_z, d_r, dp_z, dp_r, hx, hy, t, eventID)
% PLOTRADONRESULTS3D  Visualize raw vs. Radon-separated data for Z and R in 3D.
%
% Inputs:
%   d_z, dp_z : raw vs. Radon-transformed data for Z ( [Nt x Ntraces] )
%   d_r, dp_r : raw vs. Radon-transformed data for R
%   hx, hy    : x and y coordinates for each trace
%   t         : time vector (s)
%   eventID   : string identifier for the event

% Estimate amplitude for color scaling
cmax_z = 3 * rms(dp_z(:));
cmax_r = 3 * rms(dp_r(:));

figure('Name', sprintf('3D Radon Results for Event %s', eventID), ...
    'Position', [100, 100, 1200, 800], 'Color', 'w');

% Create time slice for visualization (mid-point)
time_slice = round(length(t)/2);

% Raw Z - spatial distribution at time slice
subplot(2,3,1);
scatter3(hx, hy, zeros(size(hx)), 50, d_z(time_slice,:), 'filled');
colorbar; colormap(seismic(1)); caxis([-cmax_z cmax_z]);
xlabel('X (km)'); ylabel('Y (km)'); zlabel('Amplitude');
title('Raw Z (spatial)'); set(gca, 'FontSize', 12);
view(2);

% Radon Z - spatial distribution at time slice
subplot(2,3,2);
scatter3(hx, hy, zeros(size(hx)), 50, dp_z(time_slice,:), 'filled');
colorbar; colormap(seismic(1)); caxis([-cmax_z cmax_z]);
xlabel('X (km)'); ylabel('Y (km)'); zlabel('Amplitude');
title('Radon Z (spatial)'); set(gca, 'FontSize', 12);
view(2);

% Z difference - spatial distribution at time slice
subplot(2,3,3);
scatter3(hx, hy, zeros(size(hx)), 50, d_z(time_slice,:) - dp_z(time_slice,:), 'filled');
colorbar; colormap(seismic(1)); caxis([-cmax_z cmax_z]);
xlabel('X (km)'); ylabel('Y (km)'); zlabel('Amplitude');
title('Removed Noise Z (spatial)'); set(gca, 'FontSize', 12);
view(2);

% Raw R - spatial distribution at time slice
subplot(2,3,4);
scatter3(hx, hy, zeros(size(hx)), 50, d_r(time_slice,:), 'filled');
colorbar; colormap(seismic(1)); caxis([-cmax_r cmax_r]);
xlabel('X (km)'); ylabel('Y (km)'); zlabel('Amplitude');
title('Raw R (spatial)'); set(gca, 'FontSize', 12);
view(2);

% Radon R - spatial distribution at time slice
subplot(2,3,5);
scatter3(hx, hy, zeros(size(hx)), 50, dp_r(time_slice,:), 'filled');
colorbar; colormap(seismic(1)); caxis([-cmax_r cmax_r]);
xlabel('X (km)'); ylabel('Y (km)'); zlabel('Amplitude');
title('Radon R (spatial)'); set(gca, 'FontSize', 12);
view(2);

% R difference - spatial distribution at time slice
subplot(2,3,6);
scatter3(hx, hy, zeros(size(hx)), 50, d_r(time_slice,:) - dp_r(time_slice,:), 'filled');
colorbar; colormap(seismic(1)); caxis([-cmax_r cmax_r]);
xlabel('X (km)'); ylabel('Y (km)'); zlabel('Amplitude');
title('Removed Noise R (spatial)'); set(gca, 'FontSize', 12);
view(2);
end
