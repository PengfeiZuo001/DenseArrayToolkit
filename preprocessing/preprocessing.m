function DataStruct = preprocessing(DataStruct, param)
% PREPROCESSING - Example of refactoring SNR into a separate function
%
% Usage:
%   DataStruct = preprocessing(DataStruct, param)
%
% Inputs:
%   DataStruct : the struct array with fields:
%       Waveforms.data  -> raw waveforms [Nt x 3]
%       TimeAxis.b, TimeAxis.o, TimeAxis.dt -> from SAC header
%       ...
%   param : struct of parameters (tstart, tend, sig_leader, etc.)
%
% Outputs:
%   DataStruct : updated with processed data and added fields:
%       .Waveforms.dataProcessed
%       .Waveforms.chNameProcessed
%       .TimeAxis.t_resample
%       .ProcHistory (updated logs)
%
% The function:
%   - Checks distance & keeps only [30,95] deg
%   - Finds P arrival by taupTime
%   - Calls calcSNR() to compute the signal-to-noise ratio for Z
%   - Removes DC/trend, rotates ENZ->TRZ, bandpass, chop near P, resamples
%   - Records each step in ProcHistory
%
% Author: Yunfeng Chen
% Date:   Jan. 9, 2025

%% 1. Parameter defaults
if nargin < 2, param = struct(); end
if ~isfield(param, 'tstart'),          param.tstart = -60; end
if ~isfield(param, 'tend'),            param.tend = 600;   end
if ~isfield(param, 'sig_leader'),      param.sig_leader = 30; end
if ~isfield(param, 'record_len'),      param.record_len = 120; end
if ~isfield(param, 'lows'),            param.lows = 0.02;  end
if ~isfield(param, 'highs'),           param.highs = 3;    end
if ~isfield(param, 'resample_period'), param.resample_period = 0.1; end

removeIdx = [];
disp('--- Start Preprocessing ---');
tic;

%% 2. Main loop
for n = 1:length(DataStruct)
    if mod(n, 50) == 0
        disp(['Processing trace ' num2str(n) ' / ' num2str(length(DataStruct))]);
    end
    
    %% 2.1 Basic checks: waveforms/data existence
    if ~isfield(DataStruct(n).Waveforms, 'data') || isempty(DataStruct(n).Waveforms.data)
        msgSkip = '[Warning] No waveform data found, skipping this record.';
        DataStruct(n).ProcHistory{end+1} = msgSkip;
        removeIdx(end+1) = n;
        continue;
    end

    seis = DataStruct(n).Waveforms.data;  % [Nt x 3], e.g. ENZ
    [nt,nchan] = size(seis);

    if size(seis,2) < 3
        msgSkip = '[Warning] Not enough components (<3), skipping.';
        DataStruct(n).ProcHistory{end+1} = msgSkip;
        removeIdx(end+1) = n;
        continue;
    end

    %% 2.2 Retrieve event/station info
    evla = DataStruct(n).EventInfo.evla;
    evlo = DataStruct(n).EventInfo.evlo;
    evdp = DataStruct(n).EventInfo.evdp;
    stla = DataStruct(n).StationInfo.stla;
    stlo = DataStruct(n).StationInfo.stlo;

    % distHaversine
    try
        [az, baz, degree] = distHaversine(evla, evlo, stla, stlo);
    catch ME
        msgErr = sprintf('[Error] distHaversine failed: %s', ME.message);
        DataStruct(n).ProcHistory{end+1} = msgErr;
        removeIdx(end+1) = n;
        continue;
    end
    DataStruct(n).TravelInfo.distDeg = degree;
    DataStruct(n).TravelInfo.az   = az;
    DataStruct(n).TravelInfo.baz  = baz;

    % keep only [30, 95]
    if degree < 30 || degree > 95
        msgRange = sprintf('[Info] dist=%.1f out of [30,95], removed.', degree);
        DataStruct(n).ProcHistory{end+1} = msgRange;
        removeIdx(end+1) = n;
        continue;
    end

    %% 2.3 arrival times
    try
        ptime = taupTime('prem', evdp, 'P', 'deg', degree);
        if isempty(ptime)
            msgPT = '[Warning] no P arrival found, removing.';
            DataStruct(n).ProcHistory{end+1} = msgPT;
            removeIdx(end+1) = n;
            continue;
        end
        pTime = ptime(1).time;
        DataStruct(n).TravelInfo.pTime = pTime;
        DataStruct(n).TravelInfo.rayParam  = ptime(1).rayParam / 111.1775;
    catch ME
        msgPT = sprintf('[Error] taupTime: %s', ME.message);
        DataStruct(n).ProcHistory{end+1} = msgPT;
        removeIdx(end+1) = n;
        continue;
    end

    %% 2.4 Compute SNR (move to separate function)
    % Suppose we have dt, b, o in DataStruct(n).TimeAxis
    %   dt : sampling interval
    %   b  : start time of SAC (relative)
    %   o  : event origin time in SAC reference
    if isfield(DataStruct(n).TimeAxis, 'dt'), dt = DataStruct(n).TimeAxis.dt;  else, dt = []; end
    if isfield(DataStruct(n).TimeAxis, 'b'),  b  = DataStruct(n).TimeAxis.b;   else, b  = 0;  end
    if isfield(DataStruct(n).TimeAxis, 'o'),  o  = DataStruct(n).TimeAxis.o;   else, o  = 0;  end
    if isfield(DataStruct(n).TimeAxis, 't'),  t  = DataStruct(n).TimeAxis.t;   else, t  = (0:nt-1)*dt;  end
    
    % the line below is for a specific case, wherein we request the event
    % waveform 60 sec before the P arrival,
%     p_time = 60;
    p_time = DataStruct(n).TravelInfo.pTime;
    
    try
        Z = seis(:,3);  % third column is Z
        snrVal = calcSNR(Z, dt, b, o, p_time);
        DataStruct(n).RF.snr = snrVal;
        DataStruct(n).ProcHistory{end+1} = ...
            sprintf('[Info] SNR=%.2f computed for Z.', snrVal);
    catch ME
        msgSNR = sprintf('[Warning] calcSNR failed: %s', ME.message);
        DataStruct(n).ProcHistory{end+1} = msgSNR;
        DataStruct(n).RF.snr = -999;
    end

    %% 2.5 Remove DC & trend
    seis = removeSeisDC(seis);
    seis = removeSeisTrend(seis);
    DataStruct(n).ProcHistory{end+1} = '[Info] DC & trend removed.';

    %% 2.6 Rotate ENZ -> TRZ
    try
        [seis_rot, ~] = rotateSeisENZtoTRZ(seis, baz);
        DataStruct(n).ProcHistory{end+1} = sprintf('[Info] rotateSeisENZtoTRZ done, baz=%.2f.', baz);
    catch ME
        msgRot = sprintf('[Error] rotateSeisENZtoTRZ failed: %s', ME.message);
        DataStruct(n).ProcHistory{end+1} = msgRot;
        removeIdx(end+1) = n;
        continue;
    end

    %% 2.7 Filter & cut
    seis_taper = taperSeis(seis_rot, 0.2);
    DataStruct(n).ProcHistory{end+1} = '[Info] taperSeis done.';

    if ~isempty(dt)
        seis_flt = bandpassSeis(seis_taper, dt, param.lows, param.highs, 3);
        DataStruct(n).ProcHistory{end+1} = ...
          sprintf('[Info] bandpassSeis (%.2f-%.2f Hz) done.', param.lows, param.highs);
    else
        seis_flt = seis_taper;
        DataStruct(n).ProcHistory{end+1} = '[Warning] dt not found, skip bandpassSeis.';
    end

    %% 2.8 Chop around P wave
    startCut = p_time - param.sig_leader;                % e.g. 30s before P
    finishCut= startCut + param.sig_leader + param.record_len; % total window
    [seis_cut, t_cut] = chopSeis(seis_flt, t, startCut, finishCut);
    DataStruct(n).ProcHistory{end+1} = ...
        sprintf('[Info] chopSeis around p_time=%.1f: [%.1f -> %.1f].', p_time, startCut, finishCut);

    %% 2.9 Resample
    [seis_rsp, dt_rsp, t_rsp] = resampleSeis(seis_cut, t_cut, param.resample_period);
    DataStruct(n).ProcHistory{end+1} = ...
       sprintf('[Info] resampled from dt=%.3f to %.3f.', dt, dt_rsp);

    %% 2.10 Save processed data
    DataStruct(n).Waveforms.dataProcessed    = seis_rsp;
    DataStruct(n).Waveforms.chNameProcessed  = {'T','R','Z'};
    DataStruct(n).TimeAxis.t_resample        = t_rsp;
    DataStruct(n).TimeAxis.dt_resample       = dt_rsp;
end

%% 3. Remove invalid or out-of-range records
if ~isempty(removeIdx)
    DataStruct(removeIdx) = [];
    disp(['Removed ' num2str(length(removeIdx)) ' invalid traces.']);
end

toc;
disp('--- Preprocessing completed ---');
end

function snrVal = calcSNR(Z, dt, b, o, p_time)
% calcSNR - Compute signal-to-noise ratio for the Z component
% using a simple variance ratio approach:
%   noise window: [p_time-105, p_time-5]
%   signal window: [p_time-5, p_time+5]
%
% Inputs:
%   Z      : [Nt x 1], vertical waveform
%   dt     : sampling interval
%   b, o   : from SAC (b= start time of record, o= origin time), in seconds
%   p_time : P wave arrival time (relative to origin in SAC reference)
%
% Output:
%   snrVal : ratio of var(signal) / var(noise)
%
% Author: Yunfeng Chen
% Date  : Jan. 9, 2025

if isempty(dt) || dt <= 0
    error('calcSNR:InvalidDT', 'dt must be > 0.');
end

% Convert times to sample indices:
% Index( time ) = round( ( time + (o - b) ) / dt )
% Example usage:
noiseStart = p_time - 105; 
noiseEnd   = p_time -   5;
sigStart   = p_time -   5;
sigEnd     = p_time +   5;

N = length(Z);

% Convert each to sample index
idxNoise1 = round((noiseStart + (o - b)) / dt);
idxNoise2 = round((noiseEnd   + (o - b)) / dt);
idxSig1   = round((sigStart  + (o - b)) / dt);
idxSig2   = round((sigEnd    + (o - b)) / dt);

% Bound checks
idxNoise1 = max(idxNoise1, 1); 
idxNoise2 = min(idxNoise2, N);
idxSig1   = max(idxSig1, 1);
idxSig2   = min(idxSig2, N);

if idxNoise1 >= idxNoise2 || idxSig1 >= idxSig2
    % If the windows are invalid
    error('calcSNR:InvalidWindow', ...
          'Noise or Signal window indices are out of range or overlapping.');
end

noi = Z(idxNoise1 : idxNoise2);
sig = Z(idxSig1   : idxSig2);

snrVal = var(sig) / var(noi);
end