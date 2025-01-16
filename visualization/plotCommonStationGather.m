function plotCommonStationGather(DataStruct,station)
idx = [];
for n = 1:length(DataStruct)
    if strcmp(DataStruct(n).StationInfo.sta, station)
        idx(end+1) = n;
    end
end

figure
set(gcf,'Name',['Common Station Gather: ', station],'Position',[500 500 600 1200],'Color','w');
scaling_factor = 10;
for n = idx
    % plot traces
    x = DataStruct(n).RF.itr*scaling_factor;
    t = DataStruct(n).RF.ittime;
        %     rayp = DataStruct(n).TravelInfo.rayParam;
    dist = DataStruct(n).TravelInfo.distDeg;
    p1 = plot(t,x+dist,'k','linewidth',0.1);hold on;
    % fill in the color for positive and negative phases
    thre = 0.01;
    upper = x; upper(upper<=thre) = 0;
    upper(1) = 0; upper(end) = 0;
    lower = zeros(length(x),1);
    jbfill(t(:),upper+dist,lower+dist,'r','k',1,1.0);hold on;

    upper = zeros(length(x),1);
    lower = x; lower(lower>=-thre)=0;
    lower(1) = 0; lower(end) = 0;
    jbfill(t(:),upper+dist,lower+dist,[0.17,0.17,0.17],'k',1,1.0);hold on;
    set(p1,'LineWidth',1)
    title(['Common Station Gather: ', station])
end

xlim([t(1),30])
ylim([25,100])
xlabel('Time (sec)')
ylabel('Distance (deg)')
set(gca,'fontsize',14)



