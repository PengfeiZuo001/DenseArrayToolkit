function DataStruct = fkFilter(DataStruct, gridStruct, param)
% FKFILTER  Perform FK filter-based array processing on DataStruct.
%
% Usage:
%   DataStruct = fkFilter(DataStruct, gridStruct, param)
%
% Inputs:
%   DataStruct : struct array with fields:
%       .EventInfo.evid         - Event ID (string)
%       .Waveforms.dataProcessed - Processed waveform data [Nt x 3]
%       .RF.ittime              - (not always used here, but often stored in RF)
%       .TravelInfo.distDeg     - Epicentral distance in degrees
%       .ProcHistory            - Cell array to record logs (optional)
%
%   gridStruct: struct array with grid information for station projection
%
%   param : struct with (optional) fields:
%       .lows      (default: 0.1)    - Low corner freq for bandpass (Hz)
%       .highs     (default: 1.2)    - High corner freq for bandpass (Hz)
%       .minTraces (default: 60)     - Min number of traces per event
%       .w         (default: 0.1)    - Half width of cone filter (percentage)
%       .plotFK    (default: false)  - Flag to plot FK filter results
%       .plotFKspectrum (default: false) - Flag to plot FK spectrum
%       .order     (default: postdecon) - Apply filtering before or after deconvolution
%
% Output:
%   DataStruct : The updated struct array. For each trace that belongs
%                to an event with enough traces, the field:
%       DataStruct(n).Waveforms.dataFKFiltered
%                is populated with [Nt x 3] radial/transverse (or vertical)
%                waveforms after FK-based noise removal or separation.
%
% Author:  Yunfeng Chen, Zhejiang University
% Date   : Aug. 2, 2025

%% 1. Parameter Handling & Validation
if nargin < 3 || isempty(param)
    param = struct();
end

% Define default parameters if missing
if ~isfield(param, 'lows'),       param.lows       = 0.1;   end
if ~isfield(param, 'highs'),      param.highs      = 1.2;   end
if ~isfield(param, 'minTraces'),  param.minTraces  = 60;    end
if ~isfield(param, 'w'),          param.w          = 0.1;   end
if ~isfield(param, 'plotFK'),     param.plotFK     = false; end
if ~isfield(param, 'plotFKspectrum'), param.plotFKspectrum = false; end
if ~isfield(param, 'order'),      param.order      = 'postdecon'; end

% Basic field checks (on the first element, assuming consistent struct array)
requiredTopLevel = {'EventInfo','Waveforms','TravelInfo','RF'};
for f = requiredTopLevel
    if ~isfield(DataStruct, f{1})
        error('fkFilter:MissingField',...
            'DataStruct must contain the field %s in each element.', f{1});
    end
end
if ~isfield(DataStruct(1).EventInfo, 'evid')
    error('fkFilter:MissingField',...
        'DataStruct.EventInfo must contain field: evid.');
end
if ~isfield(DataStruct(1).Waveforms, 'dataProcessed')
    error('fkFilter:MissingField',...
        'DataStruct.Waveforms must contain field: dataProcessed.');
end
if ~isfield(DataStruct(1).TravelInfo, 'distDeg')
    error('fkFilter:MissingField',...
        'DataStruct.TravelInfo must contain field: distDeg.');
end

%% 2. Identify Unique Events
eventListStruct = getEvents(DataStruct);  % retrieve unique events from the data
eventIDs = {eventListStruct.evid};

%% 3. Process Each Event
for iEvt = 1:length(eventIDs)
    eventID = eventIDs{iEvt};
    [commonEventGather, matchIndex] = getCommonEventGather(DataStruct, eventID);

    % Check if enough traces are available for FK filtering
    if length(commonEventGather) < param.minTraces
        skipMsg = sprintf('[FKFilter] Event %s skipped: only %d traces (< %d).',...
            eventID, length(commonEventGather), param.minTraces);
        commonEventGather = appendHistory(commonEventGather, skipMsg);
        DataStruct(matchIndex) = commonEventGather;
        continue;
    end

    % Display progress
    fprintf('Processing event #%d: %s with %d traces.\n', ...
        iEvt, eventID, length(commonEventGather));

    %% 3.1 Compute Distance Offsets for Each Trace
    % Use the projected station location on principle axis to define the
    % distance vector
    stationList = getStations(DataStruct);
    stlo = [stationList.stlo]';  % 台站经度
    stla = [stationList.stla]';  % 台站纬度
    [h, ~] = latlonToProjectedCoords(stlo, stla, gridStruct);
    [h, idx] = sort(h);
    commonEventGather = commonEventGather(idx);

    %% 3.2 Extract and Preprocess Waveforms
    [d_z, d_r, d_t, t, dt] = extractAndPreprocessWaveforms(commonEventGather,...
        param.lows, param.highs);
    trace_energy = rms(d_z);
    % remove traces with anomalous amplitude
    remove_idx = trace_energy > mean(trace_energy) + 2*std(trace_energy);
    d_z(:, remove_idx) = 0;
    d_r(:, remove_idx) = 0;

    switch param.order
        case 'predecon'
            %% 3.3 Handle Irregular Station Spacing for FK Filtering
            % FK filtering requires regular spacing, so interpolate to regular grid
            h_reg = linspace(min(h), max(h), length(h));
            
            % Interpolate Z component to regular grid
            d_z_reg = zeros(size(d_z));
            for i = 1:size(d_z, 1)
                d_z_reg(i, :) = interp1(h, d_z(i, :), h_reg, 'linear', 'extrap');
            end
            
            % Interpolate R component to regular grid
            d_r_reg = zeros(size(d_r));
            for i = 1:size(d_r, 1)
                d_r_reg(i, :) = interp1(h, d_r(i, :), h_reg, 'linear', 'extrap');
            end

            %% 3.4 Perform FK Filtering on Z component (regular grid)
            try
                dp_z_reg = amf_fk_dip(d_z_reg, param.w);
            catch ME
                warnMsg = sprintf('[FKFilter] Event %s, Z-component failed: %s', ...
                    eventID, ME.message);
                warning(warnMsg);
                commonEventGather = appendHistory(commonEventGather, warnMsg);
                DataStruct(matchIndex) = commonEventGather;
                continue;
            end

            %% 3.5 Perform FK Filtering on R component (regular grid)
            try
                dp_r_reg = amf_fk_dip(d_r_reg, param.w);
            catch ME
                warnMsg = sprintf('[FKFilter] Event %s, R-component failed: %s', ...
                    eventID, ME.message);
                warning(warnMsg);
                commonEventGather = appendHistory(commonEventGather, warnMsg);
                DataStruct(matchIndex) = commonEventGather;
                continue;
            end

            %% 3.6 Interpolate Results Back to Original Irregular Positions
            % Interpolate filtered Z component back to original positions
            dp_z = zeros(size(d_z));
            for i = 1:size(dp_z_reg, 1)
                dp_z(i, :) = interp1(h_reg, dp_z_reg(i, :), h, 'linear', 'extrap');
            end
            
            % Interpolate filtered R component back to original positions
            dp_r = zeros(size(d_r));
            for i = 1:size(dp_r_reg, 1)
                dp_r(i, :) = interp1(h_reg, dp_r_reg(i, :), h, 'linear', 'extrap');
            end

            %% 3.5 Store FK-Filtered Waveforms
            for n = 1:length(commonEventGather)
                % Initialize dataFKFiltered if needed
                ntData = size(d_z, 1);  % or length(t)
                commonEventGather(n).Waveforms.dataFKFiltered = zeros(ntData, 3);

                % T component (not filtered, keep original or zeros)
                commonEventGather(n).Waveforms.dataFKFiltered(:, 1) = 0;
                % R component
                commonEventGather(n).Waveforms.dataFKFiltered(:, 2) = dp_r(:, n);
                % Z component
                commonEventGather(n).Waveforms.dataFKFiltered(:, 3) = dp_z(:, n);
            end

            % Log the successful FK filtering
            successMsg = sprintf('[FKFilter] Event %s: FK filter applied (%d traces).',...
                eventID, length(commonEventGather));
            commonEventGather = appendHistory(commonEventGather, successMsg);

            % Place updated event gather back into DataStruct
            DataStruct(matchIndex) = commonEventGather;

            %% 3.6 Optional: Plot FK Filter Results and Spectrum
            if param.plotFK
                plotFKResults(d_z, d_r, dp_z, dp_r, h, t, eventID);
            end
            
            if param.plotFKspectrum
                plotFKSpectrum(d_z_reg, dp_z_reg, d_r_reg, dp_r_reg, h_reg, t, dt, eventID, param.w);
            end

        case 'postdecon'
            %% Perform FK Filtering on RF with Irregular Station Spacing Handling
            itrCell = {commonEventGather.RF};
            ittime = commonEventGather(1).RF.ittime;
            d = cell2mat(cellfun(@(rf) rf.itr, itrCell, 'UniformOutput', false));
            d(:, remove_idx) = 0;
            
            % FK filtering requires regular spacing, so interpolate to regular grid
            h_reg = linspace(min(h), max(h), length(h));
            
            % Interpolate RF data to regular grid
            d_reg = zeros(size(d));
            for i = 1:size(d, 1)
                d_reg(i, :) = interp1(h, d(i, :), h_reg, 'linear', 'extrap');
            end
            
            try
                dp_reg = amf_fk_dip(d_reg, param.w);
            catch ME
                warnMsg = sprintf('[FKFilter] Event %s, RF failed: %s', ...
                    eventID, ME.message);
                warning(warnMsg);
                commonEventGather = appendHistory(commonEventGather, warnMsg);
                DataStruct(matchIndex) = commonEventGather;
                continue;
            end

            % Interpolate filtered RF back to original irregular positions
            dp = zeros(size(d));
            for i = 1:size(dp_reg, 1)
                dp(i, :) = interp1(h_reg, dp_reg(i, :), h, 'linear', 'extrap');
            end

            for n = 1:length(commonEventGather)
                % save filtered RF
                commonEventGather(n).RF.itr = dp(:, n);
            end

            % Log the successful FK filtering
            successMsg = sprintf('[FKFilter] Event %s: FK filter applied (%d traces).',...
                eventID, length(commonEventGather));
            commonEventGather = appendHistory(commonEventGather, successMsg);

            % Place updated event gather back into DataStruct
            DataStruct(matchIndex) = commonEventGather;

            %% 3.7 Optional: Plot RF Filter Results and Spectrum for postdecon
            if param.plotFK
                plotRFResults(d, dp, h, ittime, eventID);
            end
            
            if param.plotFKspectrum
                plotRFSpectrum(d_reg, dp_reg, h_reg, ittime, dt, eventID, param.w);
            end
    end
end

end % end of fkFilter main function

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
    tmpZ = dataProc(:, 3);
    tmpR = dataProc(:, 2);
    tmpT = dataProc(:, 1);

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

function plotFKResults(d_z, d_r, dp_z, dp_r, h, t, eventID)
% PLOTFKRESULTS  Visualize raw vs. FK-filtered data for Z and R.
%
% Inputs:
%   d_z, dp_z : raw vs. FK-filtered data for Z ( [Nt x Ntraces] )
%   d_r, dp_r : raw vs. FK-filtered data for R
%   h         : distance offsets (km) for each trace
%   t         : time vector (s)
%   eventID   : string identifier for the event

% Estimate amplitude for color scaling
cmax_z = 3 * rms(dp_z(:));
cmax_r = 3 * rms(dp_r(:));

figure('Name', sprintf('FK Filter Results for Event %s', eventID), ...
    'Position', [100, 100, 1200, 800], 'Color', 'w');

% Raw Z
subplot(2, 3, 1);
imagesc(h, t, d_z);
colormap(seismic(1));
caxis([-cmax_z cmax_z]);
ylim([t(250), t(550)]);
xlabel('Distance (km)'); ylabel('Time (s)');
title('Raw Z'); set(gca, 'FontSize', 12);

% FK Filtered Z
subplot(2, 3, 2);
imagesc(h, t, dp_z);
colormap(seismic(1));
caxis([-cmax_z cmax_z]);
ylim([t(250), t(550)]);
xlabel('Distance (km)'); ylabel('Time (s)');
title('FK Filtered Z'); set(gca, 'FontSize', 12);

% Z difference
subplot(2, 3, 3);
imagesc(h, t, d_z - dp_z);
colormap(seismic(1));
caxis([-cmax_z cmax_z]);
ylim([t(250), t(550)]);
xlabel('Distance (km)'); ylabel('Time (s)');
title('Removed Noise (Z)'); set(gca, 'FontSize', 12);

% Raw R
subplot(2, 3, 4);
imagesc(h, t, d_r);
colormap(seismic(1));
caxis([-cmax_r cmax_r]);
ylim([t(250), t(550)]);
xlabel('Distance (km)'); ylabel('Time (s)');
title('Raw R'); set(gca, 'FontSize', 12);

% FK Filtered R
subplot(2, 3, 5);
imagesc(h, t, dp_r);
colormap(seismic(1));
caxis([-cmax_r cmax_r]);
ylim([t(250), t(550)]);
xlabel('Distance (km)'); ylabel('Time (s)');
title('FK Filtered R'); set(gca, 'FontSize', 12);

% R difference
subplot(2, 3, 6);
imagesc(h, t, d_r - dp_r);
colormap(seismic(1));
caxis([-cmax_r cmax_r]);
ylim([t(250), t(550)]);
xlabel('Distance (km)'); ylabel('Time (s)');
title('Removed Noise (R)'); set(gca, 'FontSize', 12);
end

function plotFKSpectrum(d_z_reg, dp_z_reg, d_r_reg, dp_r_reg, h_reg, t, dt, eventID, w)
% PLOTFKSPECTRUM  Visualize FK spectrum before and after filtering.
%
% Inputs:
%   d_z_reg, dp_z_reg : raw vs. FK-filtered data for Z on regular grid
%   d_r_reg, dp_r_reg : raw vs. FK-filtered data for R on regular grid
%   h_reg             : regular distance grid (km)
%   t                 : time vector (s)
%   dt                : time sampling (s)
%   eventID           : string identifier for the event
%   w                 : half width of cone filter (percentage)

% Compute FK spectrum for Z component
[FZ_raw, f_axis, k_axis] = computeFKSpectrum(d_z_reg, h_reg, dt);
[FZ_filt, ~, ~] = computeFKSpectrum(dp_z_reg, h_reg, dt);

% Compute FK spectrum for R component
[FR_raw, ~, ~] = computeFKSpectrum(d_r_reg, h_reg, dt);
[FR_filt, ~, ~] = computeFKSpectrum(dp_r_reg, h_reg, dt);

% Create figure
figure('Name', sprintf('FK Spectrum for Event %s', eventID), ...
    'Position', [100, 100, 1400, 1000], 'Color', 'w');

% Z component - Raw FK spectrum
subplot(2, 4, 1);
imagesc(k_axis, f_axis, log10(abs(FZ_raw)));
colormap(jet);
caxis([-4, 0]); % Adjust color scale as needed
xlabel('Wavenumber (1/km)'); ylabel('Frequency (Hz)');
title('Raw Z - FK Spectrum'); set(gca, 'FontSize', 10);
colorbar;

% Z component - Filtered FK spectrum
subplot(2, 4, 2);
imagesc(k_axis, f_axis, log10(abs(FZ_filt)));
colormap(jet);
caxis([-4, 0]);
xlabel('Wavenumber (1/km)'); ylabel('Frequency (Hz)');
title('Filtered Z - FK Spectrum'); set(gca, 'FontSize', 10);
colorbar;

% Z component - Difference (removed noise)
subplot(2, 4, 3);
FZ_diff = abs(FZ_raw) - abs(FZ_filt);
imagesc(k_axis, f_axis, FZ_diff);
colormap(jet);
caxis([-0.5, 0.5]); % Adjust color scale as needed
xlabel('Wavenumber (1/km)'); ylabel('Frequency (Hz)');
title('Z - Removed Noise (FK)'); set(gca, 'FontSize', 10);
colorbar;

% Z component - Cone filter mask overlay
subplot(2, 4, 4);
plotConeFilterMask(f_axis, k_axis, w);
xlabel('Wavenumber (1/km)'); ylabel('Frequency (Hz)');
title('Cone Filter Mask'); set(gca, 'FontSize', 10);

% R component - Raw FK spectrum
subplot(2, 4, 5);
imagesc(k_axis, f_axis, log10(abs(FR_raw)));
colormap(jet);
caxis([-4, 0]);
xlabel('Wavenumber (1/km)'); ylabel('Frequency (Hz)');
title('Raw R - FK Spectrum'); set(gca, 'FontSize', 10);
colorbar;

% R component - Filtered FK spectrum
subplot(2, 4, 6);
imagesc(k_axis, f_axis, log10(abs(FR_filt)));
colormap(jet);
caxis([-4, 0]);
xlabel('Wavenumber (1/km)'); ylabel('Frequency (Hz)');
title('Filtered R - FK Spectrum'); set(gca, 'FontSize', 10);
colorbar;

% R component - Difference (removed noise)
subplot(2, 4, 7);
FR_diff = abs(FR_raw) - abs(FR_filt);
imagesc(k_axis, f_axis, FR_diff);
colormap(jet);
caxis([-0.5, 0.5]);
xlabel('Wavenumber (1/km)'); ylabel('Frequency (Hz)');
title('R - Removed Noise (FK)'); set(gca, 'FontSize', 10);
colorbar;

% Add text information
subplot(2, 4, 8);
axis off;
text(0.1, 0.9, sprintf('Event: %s', eventID), 'FontSize', 12);
text(0.1, 0.8, sprintf('Filter width: %.1f%%', w*100), 'FontSize', 12);
text(0.1, 0.7, sprintf('Station spacing: %.1f km', mean(diff(h_reg))), 'FontSize', 12);
text(0.1, 0.6, sprintf('Frequency range: %.1f-%.1f Hz', min(f_axis), max(f_axis)), 'FontSize', 12);
text(0.1, 0.5, sprintf('Wavenumber range: %.1f-%.1f 1/km', min(k_axis), max(k_axis)), 'FontSize', 12);
end

function [F, f_axis, k_axis] = computeFKSpectrum(data, h, dt)
% COMPUTEFKSPECTRUM  Compute 2D FFT of data for FK spectrum.
%
% Inputs:
%   data : [Nt x Nx] data matrix
%   h    : distance vector (km)
%   dt   : time sampling (s)
%
% Outputs:
%   F      : 2D FFT of data (complex)
%   f_axis : frequency axis (Hz)
%   k_axis : wavenumber axis (1/km)

[Nt, Nx] = size(data);

% Compute 2D FFT
F = fft2(data);

% Shift zero frequency to center
F = fftshift(F);

% Create frequency axis
df = 1/(Nt*dt);
f_axis = (-Nt/2:Nt/2-1)*df;

% Create wavenumber axis
dh = mean(diff(h));
dk = 1/(Nx*dh);
k_axis = (-Nx/2:Nx/2-1)*dk;

% Only keep positive frequencies (since real signal spectrum is symmetric)
pos_freq_idx = f_axis >= 0;
f_axis = f_axis(pos_freq_idx);
F = F(pos_freq_idx, :);
end

function plotConeFilterMask(f_axis, k_axis, w)
% PLOTCONEFILTERMASK  Plot the cone filter mask used in FK filtering.
%
% Inputs:
%   f_axis : frequency axis (Hz)
%   k_axis : wavenumber axis (1/km)
%   w      : half width of cone filter (percentage)

Nf = length(f_axis);
Nk = length(k_axis);
mask = zeros(Nf, Nk);

% Create cone filter mask (similar to amf_fk_dip)
nw = w*Nk;
for i1 = 1:Nf
    for i2 = 1:Nk
        if i1 > (Nf/nw)*(i2-(Nk/2)) && i1 > (Nf/nw)*((Nk/2)-i2)
            mask(i1, i2) = 1;
        end
    end
end

imagesc(k_axis, f_axis, mask);
colormap(gray);
xlabel('Wavenumber (1/km)'); ylabel('Frequency (Hz)');
title('Cone Filter Mask');
end

function plotRFResults(d, dp, h, t, eventID)
% PLOTRFRESULTS  Visualize raw vs. FK-filtered receiver function data.
%
% Inputs:
%   d, dp : raw vs. FK-filtered receiver function data [Nt x Ntraces]
%   h     : distance offsets (km) for each trace
%   t     : time vector (s)
%   eventID : string identifier for the event

% Estimate amplitude for color scaling
cmax = 3 * rms(dp(:));

figure('Name', sprintf('RF FK Filter Results for Event %s', eventID), ...
    'Position', [100, 100, 1200, 600], 'Color', 'w');

% Raw RF
subplot(1, 3, 1);
imagesc(h, t, d);
colormap(seismic(1));
caxis([-cmax cmax]);
ylim([0 60])
xlabel('Distance (km)'); ylabel('Time (s)');
title('Raw RF'); set(gca, 'FontSize', 12);

% FK Filtered RF
subplot(1, 3, 2);
imagesc(h, t, dp);
colormap(seismic(1));
caxis([-cmax cmax]);
ylim([0 60])
xlabel('Distance (km)'); ylabel('Time (s)');
title('FK Filtered RF'); set(gca, 'FontSize', 12);

% RF difference (removed noise)
subplot(1, 3, 3);
imagesc(h, t, d - dp);
colormap(seismic(1));
caxis([-cmax cmax]);
ylim([0 60])
xlabel('Distance (km)'); ylabel('Time (s)');
title('Removed Noise (RF)'); set(gca, 'FontSize', 12);
end

function plotRFSpectrum(d_reg, dp_reg, h_reg, t, dt, eventID, w)
% PLOTRFSPECTRUM  Visualize receiver function data and FK spectrum before and after filtering.
%
% Inputs:
%   d_reg, dp_reg : raw vs. FK-filtered RF data on regular grid
%   h_reg         : regular distance grid (km)
%   t             : time vector (s)
%   dt            : time sampling (s)
%   eventID       : string identifier for the event
%   w             : half width of cone filter (percentage)

% Compute FK spectrum for RF data
[F_raw, f_axis, k_axis] = computeFKSpectrum(d_reg, h_reg, dt);
[F_filt, ~, ~] = computeFKSpectrum(dp_reg, h_reg, dt);

% Create figure with 2x4 layout
figure('Name', sprintf('RF FK Analysis for Event %s', eventID), ...
    'Position', [100, 100, 1600, 800], 'Color', 'w');

% Estimate amplitude for color scaling
cmax = 3 * rms(dp_reg(:));

%% Top Row: Time-Distance Domain
% 1. Raw receiver function data
subplot(2, 4, 1);
imagesc(h_reg, t, d_reg);
colormap(seismic(1));
caxis([-cmax cmax]);
ylim([0 60]);
xlabel('Distance (km)'); ylabel('Time (s)');
title('Raw RF'); set(gca, 'FontSize', 11);

% 2. Filtered receiver function data
subplot(2, 4, 2);
imagesc(h_reg, t, dp_reg);
colormap(seismic(1));
caxis([-cmax cmax]);
ylim([0 60]);
xlabel('Distance (km)'); ylabel('Time (s)');
title('FK Filtered RF'); set(gca, 'FontSize', 11);

% 3. Removed noise
subplot(2, 4, 3);
imagesc(h_reg, t, d_reg - dp_reg);
colormap(seismic(1));
caxis([-cmax cmax]);
ylim([0 60]);
xlabel('Distance (km)'); ylabel('Time (s)');
title('Removed Noise'); set(gca, 'FontSize', 11);

% 4. Cone filter mask
subplot(2, 4, 4);
plotConeFilterMask(f_axis, k_axis, w);
xlabel('Wavenumber (1/km)'); ylabel('Frequency (Hz)');
title('Cone Filter Mask'); set(gca, 'FontSize', 11);

%% Bottom Row: FK Spectrum Domain
% 5. Raw FK spectrum
subplot(2, 4, 5);
imagesc(k_axis, f_axis, log10(abs(F_raw)));
cm = flipud(colormap(jet));
colormap(cm);
caxis([0 4]);
xlabel('Wavenumber (1/km)'); ylabel('Frequency (Hz)');
title('Raw RF - FK Spectrum'); set(gca, 'FontSize', 11);
colorbar;

% 6. Filtered FK spectrum
subplot(2, 4, 6);
imagesc(k_axis, f_axis, log10(abs(F_filt)));
cm = flipud(colormap(jet));
colormap(cm);
caxis([0 4]);
xlabel('Wavenumber (1/km)'); ylabel('Frequency (Hz)');
title('Filtered RF - FK Spectrum'); set(gca, 'FontSize', 11);
colorbar;

% 7. Removed noise in FK domain
subplot(2, 4, 7);
F_diff = abs(F_raw) - abs(F_filt);
imagesc(k_axis, f_axis, F_diff);
cm = flipud(colormap(jet));
colormap(cm);
caxis([-0.5, 0.5]);
xlabel('Wavenumber (1/km)'); ylabel('Frequency (Hz)');
title('Removed Noise (FK)'); set(gca, 'FontSize', 11);
colorbar;

% 8. Key information
subplot(2, 4, 8);
axis off;
text(0.1, 0.95, sprintf('Event: %s', eventID), 'FontSize', 12, 'FontWeight', 'bold');
text(0.1, 0.85, sprintf('Filter width: %.1f%%', w*100), 'FontSize', 11);
text(0.1, 0.80, sprintf('Station spacing: %.1f km', mean(diff(h_reg))), 'FontSize', 11);
text(0.1, 0.75, sprintf('Frequency range: %.1f-%.1f Hz', min(f_axis), max(f_axis)), 'FontSize', 11);
text(0.1, 0.70, sprintf('Wavenumber range: %.1f-%.1f 1/km', min(k_axis), max(k_axis)), 'FontSize', 11);
text(0.1, 0.65, sprintf('Traces: %d', size(d_reg, 2)), 'FontSize', 11);
text(0.1, 0.60, sprintf('Time samples: %d', length(t)), 'FontSize', 11);
text(0.1, 0.55, sprintf('Time range: %.1f-%.1f s', t(1), t(end)), 'FontSize', 11);
text(0.1, 0.50, sprintf('Distance range: %.1f-%.1f km', h_reg(1), h_reg(end)), 'FontSize', 11);
text(0.1, 0.45, sprintf('Sampling rate: %.1f Hz', 1/dt), 'FontSize', 11);
text(0.1, 0.40, sprintf('Max amplitude: %.3f', max(abs(dp_reg(:)))), 'FontSize', 11);
text(0.1, 0.35, sprintf('Processing: postdecon'), 'FontSize', 11);
end
