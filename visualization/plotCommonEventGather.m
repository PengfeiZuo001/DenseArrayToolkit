function plotCommonEventGather(DataStruct, EventID, plot_type)
% PLOTSINGLEEVENTRFS - Plot receiver functions for all stations of a single event
%
% Usage:
%   plotSingleEventRFs(DataStruct, eventID, plot_type)
%
% Inputs:
%   DataStruct : struct array, each element has fields:
%       .EventInfo.evid           : event ID (string)
%       .RF.ittime (time axis)    : e.g. [Nt x 1]
%       .RF.itr (receiver func)   : e.g. [Nt x 1]
%       .TravelInfo.distDeg       : distance in degrees
%   eventID    : the event ID you want to plot
%   plot_type  : 'distance' (use distDeg as X-axis) or 'trace' (use trace index)
%
% Example:
%   plotSingleEventRFs(DataStruct, 'EV12345678', 'distance');
    %% 1) 输入检查
    if nargin < 3
        warning('plotSingleEventRFs requires 3 inputs: DataStruct, eventID, plot_type. Use trace index as default.');
        plot_type = 'trace';
    end
    %% 2) 收集满足该 eventID 的 RF 波形
    rfmatrix = [];
    distArr  = [];
    foundIndices = [];

    % 下面是完整匹配:
    isMatch = arrayfun(@(x) strcmp(x.EventInfo.evid, EventID), DataStruct);
    foundIndices = find(isMatch);

    if isempty(foundIndices)
        warning('No data found for eventID = %s', EventID);
        return;
    end

    for idx = foundIndices
        % 确保 RF.itr 存在且非空
        if isfield(DataStruct(idx), 'RF') && isfield(DataStruct(idx).RF,'itr') ...
                && ~isempty(DataStruct(idx).RF.itr)
            rfmatrix(:, end+1) = DataStruct(idx).RF.itr;          % [Nt x 1] -> add as new column
        else
            % 若这个记录没有 RF.itr, 跳过
            continue;
        end

        if isfield(DataStruct(idx), 'TravelInfo') && isfield(DataStruct(idx).TravelInfo, 'distDeg')
            distArr(end+1) = DataStruct(idx).TravelInfo.distDeg; % distance in deg
        else
            distArr(end+1) = NaN;  % 或者赋予一个默认值
        end
    end

    % 若最终 rfmatrix 为空, 说明虽然找到匹配事件, 但没有有效RF波形
    if isempty(rfmatrix)
        warning('Found %d records for eventID="%s", but none has valid RF.itr data.', ...
                 numel(foundIndices), EventID);
        return;
    end

    %% 3) 确定时间轴
    %   从第一个有 RF.ittime 的记录中获取
    firstIdx = foundIndices(1);
    if isfield(DataStruct(firstIdx), 'RF') && isfield(DataStruct(firstIdx).RF,'ittime') ...
            && ~isempty(DataStruct(firstIdx).RF.ittime)
        t = DataStruct(firstIdx).RF.ittime;
    else
        warning('No valid time axis (RF.ittime) found in the first matched record.');
        return;
    end

    [nt, nx] = size(rfmatrix);

    %% 4) 作图
    figure('Name',sprintf('RF for %s',EventID),'Color','white',...
           'Position',[200 200 1200 700]);

    switch plot_type
        case 'distance'
            wigb(rfmatrix, 2, distArr, t);
            xlabel('Distance (deg)');
        case 'trace'
            wigb(rfmatrix, 2, 1:nx, t);
            xlabel('Trace index');
    end

    ylim([t(1), 20]);    % 只显示 0~20s 区间，可按需调整
    ylabel('Time (sec)');
    set(gca, 'FontSize',14, 'LineWidth',1, 'XMinorTick','on');
    title(sprintf('Event %s : showing %d station(s)', EventID, size(rfmatrix,2)));
end
