function DataStruct = deconv(DataStruct, param)
% DECONV - Perform receiver function (RF) deconvolution, storing results in DataStruct
%          and writing logs into DataStruct(n).ProcHistory
%
% Usage:
%   DataStruct = deconv(DataStruct, param)
%
% Inputs:
%   DataStruct : struct array with fields:
%       .Waveforms.dataProcessed -> [Nt x 3] = [T, R, Z]  (from your preprocessing step)
%       .TimeAxis.t_resample     -> time axis (same length as dataProcessed)
%       .ProcHistory             -> cell array for logging (optional, but recommended)
%   param : struct with fields (default values set below):
%       .gauss       : Gauss parameter for decon (default=2.5)
%       .waterlevel  : Water-level factor (default=0.01)
%       .itmax       : Max iteration for iterative decon (default=100)
%       .minderr     : Minimum error threshold for iterative decon (default=1e-5)
%       .phaseshift  : Phase/time shift for decon (default=0)
%       .verbose     : Verbose mode for internal decon functions (default=false)
%       .radonfilter : if use filtered data with Radon transform  (default=false)
%
% Output:
%   DataStruct : the same struct array, but each valid record now has
%       DataStruct(n).RF structure that stores:
%         .wlr      - Water-level RF trace
%         .wlrms    - RMS of water-level decon residual
%         .nwl      - Additional output from makeRFwater_ammon
%         .wltime   - Time axis for water-level RF
%         .itr      - Iterative decon RF trace
%         .itrms    - Final RMS (or last iteration error)
%         .ittime   - Time axis for iterative RF
%
%   同时，为每条记录在 DataStruct(n).ProcHistory 中记录去卷积过程的日志。

%% 1) Set default parameters
if nargin < 2, param = struct(); end
if ~isfield(param, 'gauss'),      param.gauss      = 2.5;   end
if ~isfield(param, 'waterlevel'), param.waterlevel = 0.01;  end
if ~isfield(param, 'itmax'),      param.itmax      = 100;   end
if ~isfield(param, 'minderr'),    param.minderr    = 1e-5;  end
if ~isfield(param, 'phaseshift'), param.phaseshift = 5;     end
if ~isfield(param, 'verbose'),    param.verbose    = false; end
if ~isfield(param, 'radonfilter'),param.radonfilter= false; end

noresult = [];  % 用来记录无法完成反褶积的索引
tic;
disp('--- Start Deconvolution ---');

%% 2) Loop over each record in DataStruct
for n = 1:length(DataStruct)
    if mod(n,100) == 0
        fprintf('Deconvolution on trace %d/%d\n', n, length(DataStruct));
    end

    % 2.1 初始化/清空 RF 字段，避免遗留
    DataStruct(n).RF = struct();  % 若之前已有 RF 字段，可先清空后再写

    % 2.2 一些存在性检查
    if ~isfield(DataStruct(n).Waveforms, 'dataProcessed') || ...
            isempty(DataStruct(n).Waveforms.dataProcessed)
        warnMsg = sprintf('[Deconv] No processed waveforms for trace %d, skip.', n);
        warning(warnMsg);
        DataStruct=appendHistory(DataStruct, n, warnMsg);
        noresult(end+1) = n;
        continue;
    else
        seisPRZ = DataStruct(n).Waveforms.dataProcessed;  % e.g., [Nt x 3] => [T,R,Z]
        if size(seisPRZ,2) < 3
            warnMsg = sprintf('[Deconv] dataProcessed has <3 comps for trace %d, skip.', n);
            warning(warnMsg);
            DataStruct=appendHistory(DataStruct, n, warnMsg);
            noresult(end+1) = n;
            continue;
        end
    end

    if param.radonfilter
        if ~isfield(DataStruct(n).Waveforms, 'dataRadonFiltered') || ...
                isempty(DataStruct(n).Waveforms.dataRadonFiltered)
            warnMsg = sprintf('[Deconv] No Radon filtered waveforms for trace %d, skip.', n);
            warning(warnMsg);
            DataStruct=appendHistory(DataStruct, n, warnMsg);
            noresult(end+1) = n;
            continue;
        else

            seisPRZ = DataStruct(n).Waveforms.dataRadonFiltered;  % e.g., [Nt x 3] => [T,R,Z]
            if size(seisPRZ,2) < 3
                warnMsg = sprintf('[Deconv] dataRadonFiltered has <3 comps for trace %d, skip.', n);
                warning(warnMsg);
                DataStruct=appendHistory(DataStruct, n, warnMsg);
                noresult(end+1) = n;
                continue;
            end
        end
    end


    % 检查时间轴
    if ~isfield(DataStruct(n).TimeAxis, 't_resample') || ...
            isempty(DataStruct(n).TimeAxis.t_resample)
        warnMsg = sprintf('[Deconv] TimeAxis.t_resample missing/empty for trace %d, skip.', n);
        warning(warnMsg);
        DataStruct=appendHistory(DataStruct, n, warnMsg);
        noresult(end+1) = n;
        continue;
    end
    t = DataStruct(n).TimeAxis.t_resample;
    if length(t) ~= size(seisPRZ,1)
        warnMsg = sprintf('[Deconv] length(t)=%d but size(seisPRZ,1)=%d mismatch, skip trace %d.', ...
            length(t), size(seisPRZ,1), n);
        warning(warnMsg);
        DataStruct=appendHistory(DataStruct, n, warnMsg);
        noresult(end+1) = n;
        continue;
    end

    % 2.3 准备 R, Z 分量
    R = seisPRZ(:,2);  % 假设: T=1, R=2, Z=3
    Z = seisPRZ(:,3);
    dt = t(2) - t(1);
    nt = length(t);

    %% 2.4 Perform water-level decon
    try
        [wlr, wlrms, nwl] = makeRFwater_ammon( ...
            R, Z, param.phaseshift, dt, nt, param.waterlevel, param.gauss, param.verbose);

        % 写入 DataStruct(n).RF
        DataStruct(n).RF.wlr   = wlr(:);
        DataStruct(n).RF.wlrms = wlrms;
        DataStruct(n).RF.nwl   = nwl;
        DataStruct(n).RF.wltime= (dt*(0:nt-1) - param.phaseshift)';

        % 记录日志
        logMsg = sprintf('[Deconv] Water-level decon success (gauss=%.2f, wlevel=%.3f)', ...
            param.gauss, param.waterlevel);
        DataStruct=appendHistory(DataStruct, n, logMsg);

    catch ME
        warnMsg = sprintf('[Deconv] Trace %d: makeRFwater_ammon failed: %s', n, ME.message);
        warning(warnMsg);
        DataStruct=appendHistory(DataStruct, n, warnMsg);
        noresult(end+1) = n;
        continue;
    end

    %% 2.5 Perform iterative decon
    try
        [itr, itrms] = makeRFitdecon_la_norm( ...
            R, Z, dt, nt, param.phaseshift, param.gauss, param.itmax, param.minderr);
    catch ME
        warnMsg = sprintf('[Deconv] Trace %d: iterative decon failed: %s', n, ME.message);
        warning(warnMsg);
        DataStruct=appendHistory(DataStruct, n, warnMsg);
        noresult(end+1) = n;
        continue;
    end

    if isempty(itrms)
        warnMsg = sprintf('[Deconv] Trace %d: empty itrms => no iterative RF.', n);
        warning(warnMsg);
        DataStruct=appendHistory(DataStruct, n, warnMsg);
        noresult(end+1) = n;
    else
        DataStruct(n).RF.itr   = itr(:);
        DataStruct(n).RF.itrms = itrms(end);
        DataStruct(n).RF.ittime= (dt*(0:nt-1) - param.phaseshift)';

        % 记录日志
        logMsg = sprintf('[Deconv] Iter decon success (gauss=%.2f, itmax=%d, minderr=%.1e)', ...
            param.gauss, param.itmax, param.minderr);
        DataStruct=appendHistory(DataStruct, n, logMsg);
    end
end

%% 3) 处理无效记录（可选）
if ~isempty(noresult)
    % 若你想在 MVP 阶段直接删除这些无结果记录，可取消注释:
    DataStruct(noresult) = [];
    fprintf('Removed %d traces with no valid decon result.\n', length(noresult));
end

toc;
disp('--- Deconvolution completed ---');
end


%% 辅助函数：appendHistory
function DataStruct=appendHistory(DataStruct, idx, msg)
% APPENDHISTORY - a small helper to append a message to DataStruct(idx).ProcHistory
% If DataStruct(idx).ProcHistory doesn't exist, we create it as a cell array.
if ~isfield(DataStruct(idx), 'ProcHistory') || isempty(DataStruct(idx).ProcHistory)
    DataStruct(idx).ProcHistory = {msg};
else
    DataStruct(idx).ProcHistory{end+1} = msg;
end
end