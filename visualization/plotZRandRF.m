function plotZRandRF(commonEventGather)
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
    tmpRF = commonEventGather(iTr).RF.itr;

%     % Taper (5 sec at start/end) - adapt function arguments as needed
%     tmpZ = taper(tmpZ, 5, 5);
%     tmpR = taper(tmpR, 5, 5);
%     tmpT = taper(tmpT, 5, 5);
% 
%     % Bandpass filter
%     tmpZ = bandpassSeis(tmpZ, dt, lows, highs, 3);
%     tmpR = bandpassSeis(tmpR, dt, lows, highs, 3);
%     tmpT = bandpassSeis(tmpT, dt, lows, highs, 3);


    % Store in output arrays
    d_z(:, iTr) = tmpZ;
    d_r(:, iTr) = tmpR;
    d_t(:, iTr) = tmpT;
    d_rf(:,iTr) = tmpRF;
end
h = 1:size(d_rf,2);
t = commonEventGather(1).TimeAxis.t;
ittime=commonEventGather(1).RF.ittime;

cmax_z = 1 * rms(d_z(:));
cmax_r = 1 * rms(d_r(:));
cmax_rf =1 * rms(d_rf(:)); 
figure
set(gcf,'Position',[0 0 1500 400],'Color','w')
% Raw Zs
subplot(1,3,1);
imagesc(h, t, d_z);
colormap(seismic(1));
caxis([-cmax_z cmax_z]);
ylim([t(1)+120,t(1)+150])
% wigb(d_z,1,h,t)
xlabel('Trace #'); ylabel('Time (s)');
title('Raw Veritcal Component'); set(gca, 'FontSize', 12);

% Raw Rs
subplot(1,3,2);
imagesc(h, t, d_r);
colormap(seismic(1));
caxis([-cmax_z cmax_z]);
ylim([t(1)+120,t(1)+150])
% wigb(d_z,1,h,t)
xlabel('Trace #'); ylabel('Time (s)');
title('Raw Radial Component'); set(gca, 'FontSize', 12);

% Raw RFs

subplot(1,3,3);
imagesc(h, ittime, d_rf);
colormap(seismic(1));
caxis([-cmax_rf cmax_rf]);
ylim([-5 30])
% wigb(d_z,1,h,t)
xlabel('Trace #'); ylabel('Time (s)');
title('Raw Receiver Functions'); set(gca, 'FontSize', 12);