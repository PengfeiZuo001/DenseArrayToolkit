function plotCCPResults(ccpResult,slon,slat)


RayMatrix_all = ccpResult.rayMatrix;
F = ccpResult.rf;

%% plot raypath

figure
% set(gcf,'Position',[100 100 1200 800])
set(gcf,'Position',[0 0 1600 800],'Color','w')
for n=1:5:size(RayMatrix_all,2)
    lattmp=RayMatrix_all(:,n,3);
    lontmp=RayMatrix_all(:,n,4);
    ztmp=RayMatrix_all(:,n,2);
    atmp=RayMatrix_all(:,n,1);
%     scatter(lattmp,ztmp,10,atmp,'filled'); hold on;
    scatter3(lontmp,lattmp,ztmp,30,atmp,'filled'); hold on;
end
% viscircles(centers,radii); hold on;
% plot station location
scatter(slon,slat,200,'b^','filled')
set(gca,'Zdir','reverse')
xlabel(['Longitude', char(176)]);
ylabel(['Latitude', char(176)]);
zlabel('Depth (km)');
set(gca,'Xdir','reverse')
set(gca,'Ydir','reverse')

set(gca,'FontSize',16)
colormap(seismic(3))
caxis([-0.5 0.5])
colorbar
% xlim([91 95])
% ylim([37 39])
zlim([0 800])
view(176.0732,52.4066)


end