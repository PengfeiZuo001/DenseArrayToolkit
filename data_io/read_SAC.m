function DataStruct = read_SAC(dataFolder, maxFiles)
% function to read in sac data and save as a struct

if nargin < 2
    maxFiles = Inf;  % 默认读取所有
end

sacfiles = dir(fullfile(dataFolder, '*Z.*.SAC'));
nFiles = min(length(sacfiles), maxFiles);
DataStruct = [];

for i = 1:nFiles
    % 显示简单进度
    if mod(i, round(nFiles/10)) == 0
        disp(['Reading file ' num2str(i) ' / ' num2str(nFiles)]);
    end

    % 初始化结构
    tmpStruct = initEmptyStruct();

    %---------- 1) 读取 Z 分量 ----------
    Zfile = fullfile(sacfiles(i).folder, sacfiles(i).name);
    [tZ, dataZ, hdrZ] = fget_sac(Zfile);

    % 如果 Z 文件都读取失败，可在此做判断 (此处假设 fget_sac 无异常)
    %----------

    %---------- 2) 生成对应的 N/E 文件名 ----------
    Zcomp = deblank(char(hdrZ.stations.kcmpnm)); % SAC 中的通道名
    Ncomp = strrep(Zcomp,'Z','N');
    Ecomp = strrep(Zcomp,'Z','E');

    Nfile = strrep(Zfile, Zcomp, Ncomp);
    Efile = strrep(Zfile, Zcomp, Ecomp);

    % 检查 N/E 文件是否存在
    missingComponents = false;
    missingMsg = '';
    if ~isfile(Nfile) || ~isfile(Efile)
        missingComponents = true;
        missingMsg = sprintf('[Warning] Missing N or E file for Zfile = %s. Skipped.\n', Zfile);
        warning(missingMsg);  % 也以 warning 形式输出
        % 在ProcHistory记录
        tmpStruct.ProcHistory{end+1} = missingMsg; 
        % 视需求可选择 continue 或仅保留Z分量
        continue  
    end

    %---------- 3) 读取 N/E 分量 ----------
    [tN, dataN, hdrN] = fget_sac(Nfile);
    [tE, dataE, hdrE] = fget_sac(Efile);

    %---------- 4) 计算三分量的绝对起始时间 (UTC) ----------
    % 建议统一以Z分量为主，若 hdrN, hdrE 与 hdrZ 不一致，可进一步核对。
    DateNumberZ = sacHeaderToDateNum(hdrZ);  % 你需要自己实现的辅助函数 
    % 它把 nzyear, nzjday, nzhour, nzmin, nzsec, nzmsec 等转换成 MATLAB datenum
    % 或者可直接在此写类似之前的 jday_to_day + datenum 逻辑

    % 绝对起始/结束时间(UTC): 
    %   (Z分量) = (绝对参考时刻) + (hdrZ.times.b / 86400)

    tmpStruct.TimeAxis.dt = hdrZ.times.delta; 
    tmpStruct.TimeAxis.b = hdrZ.times.b; 
    tmpStruct.TimeAxis.e = hdrZ.times.e; 
    tmpStruct.TimeAxis.o = hdrZ.times.o; 
    % 需要至少记录o或者P time中的一个，否则无法恢复绝对时间轴
    if tmpStruct.TimeAxis.o == -12345
        tmpStruct.TimeAxis.o = 0;
    end
    tmpStruct.TimeAxis.startTimeUTC = DateNumberZ + hdrZ.times.b/86400; 
    tmpStruct.TimeAxis.endTimeUTC   = DateNumberZ + hdrZ.times.e/86400; 
    tmpStruct.EventInfo.orginTimeUTC = DateNumberZ + tmpStruct.TimeAxis.o/86400;

    % 绝对时间向量
    tmpStruct.TimeAxis.t = (0:length(dataZ)-1)' * hdrZ.times.delta + tmpStruct.TimeAxis.b;

    %---------- 5) 检查三分量长度并截断到最短 ----------
    lenZ = length(dataZ);
    lenN = length(dataN);
    lenE = length(dataE);
    minLen = min([lenZ, lenN, lenE]);

    if (lenZ ~= lenN) || (lenN ~= lenE)
        truncateMsg = sprintf(['[Info] Station %s has length mismatch among Z,N,E. ', ...
                               'Truncating all to minLen = %d.\n'], ...
                              deblank(char(hdrZ.station.kstnm)), minLen);
        warning(truncateMsg);  % 也可根据需要用warning或disp
        tmpStruct.ProcHistory{end+1} = truncateMsg;
    end

    dataZ = dataZ(1:minLen);
    dataN = dataN(1:minLen);
    dataE = dataE(1:minLen);

    % 保证时间轴t也截断到 minLen
    tmpStruct.TimeAxis.t = tmpStruct.TimeAxis.t(1:minLen);

    % 若有必要，也可重新计算 endTimeUTC:
    %  endTimeUTC = startTimeUTC + (minLen-1)*dt/86400
    tmpStruct.TimeAxis.endTimeUTC = tmpStruct.TimeAxis.startTimeUTC + ...
                                    (minLen-1)*hdrZ.times.delta / 86400;

    %---------- 6) 组合三分量数据 ----------
    tmpStruct.Waveforms.data = [dataE, dataN, dataZ];  % [minLen x 3]
    tmpStruct.Waveforms.chName = {Ecomp, Ncomp, Zcomp};  % 记录通道名（可选）

    %---------- 7) 补充事件与台站信息 ----------
    % 下面从 hdrZ.station / hdrZ.event 中读数据
    % 若多分量头信息不一致，可做进一步核对
    tmpStruct.StationInfo.sta    = deblank(char(hdrZ.station.kstnm));
    tmpStruct.StationInfo.stla    = hdrZ.station.stla;
    tmpStruct.StationInfo.stlo    = hdrZ.station.stlo;
    tmpStruct.StationInfo.stel   = hdrZ.station.stel;
    tmpStruct.StationInfo.network= deblank(char(hdrZ.stations.knetwk));

    tmpStruct.EventInfo.evla = hdrZ.event.evla;
    tmpStruct.EventInfo.evlo = hdrZ.event.evlo;
    tmpStruct.EventInfo.evdp = hdrZ.event.evdp;
    tmpStruct.EventInfo.evid = datestr(tmpStruct.EventInfo.orginTimeUTC,'yyyymmddHHMMSS');
    
    % 如果还需年、月、日等可继续写

    %---------- 8) 记录文件路径 ----------
    tmpStruct.Header.filenameZ = Zfile;
    tmpStruct.Header.filenameN = Nfile;
    tmpStruct.Header.filenameE = Efile;

    %---------- 9) ProcHistory 初始化 ----------
    initialLog = sprintf('Data read from SAC (E,N,Z) on %s', datestr(now, 31));
    tmpStruct.ProcHistory{end+1} = initialLog;

    % 若想检查多分量头信息差异，示例:
    [inconsistency, detailMsg] = checkHeaderInconsistency(hdrZ, hdrN, hdrE);
    if inconsistency
       warning(detailMsg);
       tmpStruct.ProcHistory{end+1} = detailMsg;
       % 根据需求也可在这里决定是否 continue
    end

    %---------- 10) 存入输出结构 ----------
    DataStruct = [DataStruct, tmpStruct];
end
end

%% 初始化空结构
function tmpStruct = initEmptyStruct()
tmpStruct = struct( ...
        'Waveforms',   [], ...  % 波形数据: Nt x Nch 或 cell 数组
        'TimeAxis',    [], ...  % 绝对时间、采样率等
        'StationInfo', [], ...  % 台站信息
        'EventInfo',   [], ...  % 震源信息
        'Header',      [], ...  % 文件路径/头信息
        'RF',          [], ...  % 接收函数
        'TravelInfo',  [], ...  % 这里存 dist, az, baz, pTime, rayParam 等
        'ProcHistory', {{}} ... % 用 cell 记录处理历史
    );
end
%% 示例: 把 SAC header -> datenum
function dateNum = sacHeaderToDateNum(hdr)
% 由 hdr 里的 nzyear, nzjday, nzhour, nzmin, nzsec, nzmsec 推算绝对UTC datenum
% 注意: 需要你自己确保 hdr 结构中字段名对应
    nzyear = hdr.event.nzyear;
    nzjday = hdr.event.nzjday;
    nzhour = hdr.event.nzhour;
    nzmin  = hdr.event.nzmin;
    nzsec  = hdr.event.nzsec;
    nzmsec = hdr.event.nzmsec;

    if any([nzyear, nzjday, nzhour, nzmin, nzsec] == -12345)
        % -12345 通常表示 SAC 头信息无效
        % 记录到ProcHistory可在外部实现，或这里直接抛错
        error('Invalid SAC header: missing time info.');
    end

    [nzday, nzmonth] = jday_to_day(nzyear, nzjday);  % 你已有的函数
    date_str = sprintf('%04d-%02d-%02d %02d:%02d:%02d',...
                       nzyear, nzmonth, nzday, nzhour, nzmin, nzsec);

    dateNum = datenum(date_str, 'yyyy-mm-dd HH:MM:SS');
    % 加上毫秒
    dateNum = dateNum + (nzmsec / 1000) / 86400;
end

%% 示例: 检查多分量头信息是否不一致
function [inconsistency, detailMsg] = checkHeaderInconsistency(hdrZ, hdrN, hdrE)
% 只是示例: 可比较台站坐标、事件坐标、Delta等是否不同
inconsistency = false;
detailMsg = '';

% 比如：若 station.stla 不同，就视为不一致
if abs(hdrZ.station.stla - hdrN.station.stla) > 1e-5 ...
        || abs(hdrZ.station.stla - hdrE.station.stla) > 1e-5
    inconsistency = true;
    detailMsg = sprintf('[Warning] Station lat mismatch among Z,N,E files: %f vs %f vs %f\n',...
                        hdrZ.station.stla, hdrN.station.stla, hdrE.station.stla);
end

% 如果还想比较 station.stlo, event.evla, times.delta 等等，都可在此添加
% ...
end
