function DataStruct = radonTransform(DataStruct, param)
% RADONTRANSFORM - Perform Radon Transform-based array processing on DataStruct
%
% Usage:
%   DataStruct = radonTransform(DataStruct, param)
%
% Inputs:
%   DataStruct - struct array with fields:
%       .EventInfo.evid        - Event ID (string)
%       .Waveforms.dataProcessed - Processed waveform data [Nt x 3]
%       .RF.ittime            - Time axis for iterative RF [Nt x 1]
%       .TravelInfo.distDeg    - Distance in degrees [scalar]
%   param - struct with optional fields:
%       .lows     - Low corner frequency (Hz) (default: 0.1)
%       .highs    - High corner frequency (Hz) (default: 1.2)
%       .pmax     - Maximum slowness (s/km) (default: 0.05)
%       .pmin     - Minimum slowness (s/km) (default: -0.05)
%       .minTraces - Minimum number of traces per event (default: 60)
%       .N1       - CG iterations (default: 30)
%       .N2       - Update weight for sparse solution (default: 1)
%       .plotRadon - Whether to plot Radon results (default: false)
%
% Outputs:
%   DataStruct - Updated struct array with Radon-transformed waveforms
%
% Author: Yunfeng Chen
% Date: Jan. 12, 2025

    %% 1. Input Parameter Handling and Validation
    if nargin < 2 || isempty(param)
        param = struct();
    end
    % Set default parameters
    if ~isfield(param, 'lows'),      param.lows      = 0.1;    end
    if ~isfield(param, 'highs'),     param.highs     = 1.2;    end
    if ~isfield(param, 'pmax'),      param.pmax      = 0.05;   end
    if ~isfield(param, 'pmin'),      param.pmin      = -0.05;  end
    if ~isfield(param, 'minTraces'), param.minTraces = 60;   end
    if ~isfield(param, 'N1'),        param.N1        = 30;     end
    if ~isfield(param, 'N2'),        param.N2        = 1;      end
    if ~isfield(param, 'plotRadon'), param.plotRadon = false; end

    % Validate DataStruct
    requiredFields = {'EventInfo', 'Waveforms', 'TravelInfo', 'RF'};
    for fld = requiredFields
        if ~isfield(DataStruct, fld{1})
            error('radonTransform:MissingField', 'DataStruct must contain the field: %s.', fld{1});
        end
    end
    if ~isfield(DataStruct(1).EventInfo, 'evid')
        error('radonTransform:MissingField', 'DataStruct.EventInfo must contain the field: evid.');
    end
    if ~isfield(DataStruct(1).Waveforms, 'dataProcessed')
        error('radonTransform:MissingField', 'DataStruct.Waveforms must contain the field: dataProcessed.');
    end
    if ~isfield(DataStruct(1).TravelInfo, 'distDeg')
        error('radonTransform:MissingField', 'DataStruct.TravelInfo must contain the field: distDeg.');
    end

    %% 2. Prepare Parameters
    p = linspace(param.pmin, param.pmax, 20);  % Slowness vector
    np = length(p);
    Param.p = p;
    Param.N1 = param.N1;
    Param.N2 = param.N2;

    %% 3. Extract Unique Events
    eventList = getEvents(DataStruct);

    %% 4. Process Each Event
    for i = 1:length(eventList)
        eventID = eventList{i};
        [CommonEventGather, matchIndex] = getCommonEventGather(DataStruct, eventID);
        
        % Criteria: at least minTraces traces
        if length(CommonEventGather) < param.minTraces
            CommonEventGather=appendHistory(CommonEventGather, sprintf('[RadonTransform] Event %s skipped: only %d traces.', eventID, length(CommonEventGather)));
            DataStruct(matchIndex) = CommonEventGather;
            continue;
        end

        % Display progress
        fprintf('Processing event# %d: %s with %d traces.\n', i, eventID, length(CommonEventGather));

        %% 4.1 Extract and Preprocess Waveforms
        [d_z, d_r, d_t, t, dt] = extractAndPreprocessWaveforms(CommonEventGather, param.lows, param.highs);

        %% 4.2 Compute Distance and Setup Radon Parameters
        dist = [];
        for n = 1:length(CommonEventGather)
            dist(end+1) = CommonEventGather(n).TravelInfo.distDeg*(2*pi*6371/360);
        end

        h = dist - min(dist);
        Param.h = h;
        Param.v = 1./p;  % Slowness inversion
        Param.nt = length(t);
        Param.dt = dt;
        Param.type = 1;
        ma = zeros(Param.nt, np);  % Preallocate for Radon transform

        %% 4.3 Perform Radon Transform on Z Component
        try
            mi_z = yc_pcg(@radon_op, Param, d_z, zeros(size(ma)), param.N1, param.N2, 1);
            dp_z = radon_op(mi_z, Param, 1);
        catch ME
            warnMsg = sprintf('[RadonTransform] Event %s: Radon transform (Z) failed: %s', eventID, ME.message);
            warning(warnMsg);
            CommonEventGather=appendHistory(CommonEventGather, sprintf('[RadonTransform] Event %s skipped: only %d traces.', eventID, length(CommonEventGather)));
            DataStruct(matchIndex) = CommonEventGather;
            continue;
        end

        %% 4.4 Perform Radon Transform on R Component
        try
            mi_r = yc_pcg(@radon_op, Param, d_r, zeros(size(ma)), param.N1, param.N2, 1);
            dp_r = radon_op(mi_r, Param, 1);
        catch ME
            warnMsg = sprintf('[RadonTransform] Event %s: Radon transform (R) failed: %s', eventID, ME.message);
            warning(warnMsg);
            appendHistory(CommonEventGather, warnMsg)
            DataStruct(matchIndex) = CommonEventGather;
            continue;
        end
        %% 4.5 Perform Radon Transform on T Component
%         try
%             mi_t = yc_pcg(@radon_op, Param, d_t, zeros(size(ma)), param.N1, param.N2, 1);
%             dp_t = radon_op(mi_t, Param, 1);
%         catch ME
%             warnMsg = sprintf('[RadonTransform] Event %s: Radon transform (T) failed: %s', eventID, ME.message);
%             warning(warnMsg);
%             appendHistory(CommonEventGather, warnMsg)
%             DataStruct(matchIndex) = CommonEventGather;
%             continue;
%         end
        %% 4.6 Save Radon-Processed Data Back to DataStruct
        for n = 1:length(CommonEventGather)
%             CommonEventGather(n).Waveforms.dataRadonFiltered(:,1) = dp_t(:,n);
            CommonEventGather(n).Waveforms.dataRadonFiltered(:,1) = zeros(size(dp_z,1),1);
            CommonEventGather(n).Waveforms.dataRadonFiltered(:,2) = dp_r(:,n);
            CommonEventGather(n).Waveforms.dataRadonFiltered(:,3) = dp_z(:,n);
        end
        CommonEventGather=appendHistory(CommonEventGather, sprintf('[RadonTransform] Radon transform applied to trace %d.', matchIndex(n)));
        DataStruct(matchIndex) = CommonEventGather;

        %% 4.6 Optional: Plot Radon Results
        if param.plotRadon
            plotRadonResults(d_z, d_r, dp_z, dp_r, h, t, eventID);
        end
    end
    end

    %% ----------------------- Subfunctions -----------------------
    function [d_z, d_r, d_t, t, dt] = extractAndPreprocessWaveforms(CommonEventGather, lows, highs)
    % EXTRACTANDPREPROCESSWAVEFORMS - Extract and preprocess Z and R components
        numTraces = length(CommonEventGather);
        dataSample = CommonEventGather(1).Waveforms.dataProcessed;
        [nt, ~] = size(dataSample);
        d_z = zeros(nt, numTraces);
        d_r = zeros(nt, numTraces);
        d_t = zeros(nt, numTraces);
        dt = CommonEventGather(1).TimeAxis.t_resample(2) - CommonEventGather(1).TimeAxis.t_resample(1);
        for n = 1:numTraces
            dataProcessed = CommonEventGather(n).Waveforms.dataProcessed;
            d_z(:,n) = taper(dataProcessed(:,3), 5, 5);
            d_z(:,n) = bandpassSeis(d_z(:,n), dt, lows, highs, 3);
            
            d_r(:,n) = taper(dataProcessed(:,2), 5, 5);
            d_r(:,n) = bandpassSeis(d_r(:,n), dt, lows, highs, 3);

            d_t(:,n) = taper(dataProcessed(:,1), 5, 5);
            d_t(:,n) = bandpassSeis(d_t(:,n), dt, lows, highs, 3);
        end
        
        t = CommonEventGather(1).TimeAxis.t_resample;
        dt = t(2) - t(1);
    end

    function CommonEventGather = appendHistory(CommonEventGather, msg)
    % APPENDHISTORY - Append a message to ProcHistory of DataStruct at index idx
    for i = 1:length(CommonEventGather)
        if isempty(CommonEventGather(i).ProcHistory)
            
        else
            CommonEventGather(i).ProcHistory{end+1} = msg;
        end
    end
    end

    function plotRadonResults(d_z, d_r, dp_z, dp_r, h, t, eventID)
    % PLOTRADONRESULTS - Plot raw and Radon-transformed results for a given event
        % Define color scaling based on RMS of Radon data
        cmax_z = 3 * rms(dp_z(:));
        cmax_r = 3 * rms(dp_r(:));
        
        figure('Name', sprintf('Radon Results for Event %s', eventID), 'Position', [100 100 1200 800], 'Color', 'w');
        
        subplot(2,3,1);
        imagesc(h, t, d_z); 
        colormap(seismic(1)); 
        caxis([-cmax_z cmax_z]);
        xlabel('Distance (km)'); ylabel('Time (sec)');
        title('Raw Z'); set(gca, 'FontSize', 12);
        
        subplot(2,3,2);
        imagesc(h, t, dp_z); 
        colormap(seismic(1)); 
        caxis([-cmax_z cmax_z]);
        xlabel('Distance (km)'); ylabel('Time (sec)');
        title('Radon Z'); set(gca, 'FontSize', 12);
        
        subplot(2,3,3);
        imagesc(h, t, d_z - dp_z); 
        colormap(seismic(1)); 
        caxis([-cmax_z cmax_z]);
        xlabel('Distance (km)'); ylabel('Time (sec)');
        title('Removed Noise (Z)'); set(gca, 'FontSize', 12);
        
        subplot(2,3,4);
        imagesc(h, t, d_r); 
        colormap(seismic(1)); 
        caxis([-cmax_r cmax_r]);
        xlabel('Distance (km)'); ylabel('Time (sec)');
        title('Raw R'); set(gca, 'FontSize', 12);
        
        subplot(2,3,5);
        imagesc(h, t, dp_r); 
        colormap(seismic(1)); 
        caxis([-cmax_r cmax_r]);
        xlabel('Distance (km)'); ylabel('Time (sec)');
        title('Radon R'); set(gca, 'FontSize', 12);
        
        subplot(2,3,6);
        imagesc(h, t, d_r - dp_r); 
        colormap(seismic(1)); 
        caxis([-cmax_r cmax_r]);
        xlabel('Distance (km)'); ylabel('Time (sec)');
        title('Removed Noise (R)'); set(gca, 'FontSize', 12);
    end
