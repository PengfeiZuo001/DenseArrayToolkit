function plotWaveforms(DataStruct,trace_index)
% plot three components seismogram and RFs
t = DataStruct(trace_index).TimeAxis.t_resample;
T = DataStruct(trace_index).Waveforms.dataProcessed(:,1);
R = DataStruct(trace_index).Waveforms.dataProcessed(:,2);
Z = DataStruct(trace_index).Waveforms.dataProcessed(:,3);
ittime = DataStruct(trace_index).RF.ittime;
itr = DataStruct(trace_index).RF.itr;
wltime = DataStruct(trace_index).RF.wltime;
wlr = DataStruct(trace_index).RF.wlr;

% Set the current data value
figure;
set(gcf,'Position',[10 10 1000 800],'Color','w');
subplot(511)
plot(t,R);
xlim([t(1) t(1) + 120])
% set(gca,'Visible','off')
subplot(512)
plot(t,T);
xlim([t(1) t(1) + 120])
subplot(513)
plot(t,Z);
xlim([t(1) t(1) + 120])

subplot(514)
x = wlr;
t = wltime;
thre = 0.01;
upper = x; upper(upper<=thre) = 0;
upper(1) = 0; upper(end) = 0;
lower = zeros(length(x),1);
plot(t,x);hold on;
jbfill(t,upper,lower,'r','k',1,1.0);hold off;
xlim([t(1) 30])
ylim([-0.25 0.5])
xlabel('Time (sec)')

subplot(515)
x = itr;
t = ittime;
thre = 0.01;
upper = x; upper(upper<=thre) = 0;
upper(1) = 0; upper(end) = 0;
lower = zeros(length(x),1);
plot(t,x);hold on;
jbfill(t,upper,lower,'r','k',1,1.0);hold off;
xlim([t(1) 30])
ylim([-0.25 0.5])
xlabel('Time (sec)')
