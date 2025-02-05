function plotCCPResults(ccpResult,slon,slat)


RayMatrix_all = ccpResult.rayMatrix;
F = ccpResult.rf;

%% plot raypath

figure
% set(gcf,'Position',[100 100 1200 800])
set(gcf,'Position',[0 0 1600 800],'Color','w')
for n=1:50:size(RayMatrix_all,2)
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
% set(gca,'Xdir','reverse')
% set(gca,'Ydir','reverse')

set(gca,'FontSize',16)
colormap(seismic(3))
caxis([-0.5 0.5])
colorbar
% xlim([91 95])
% ylim([37 39])
zlim([0 120])
view([140,20])

%% plot ccp results
% need to adjust 
lon1 = 95.29;
lat1 = 41.51;
lon2 = 96.18;
lat2 = 39.49;
dz = 0.5;
zmax = 800;
dep0 = 0:dz:zmax;

nlatlon = 100;
[latp,lonp] = gcwaypts(lat1,lon1,lat2,lon2,nlatlon);
[deg0,az0]= distance(lat1,lon1,latp,lonp);
% degree to distance
dist0 = deg0*2*pi*6371/360;

% plot the geometry of RFs

LAT_profile1 = repmat(latp',length(dep0),1);
LON_profile1 = repmat(lonp',length(dep0),1);
DEP_profile1 = repmat(dep0',1,length(latp));
DIST_profile1 = repmat(dist0',length(dep0),1);
Vprofile1 = F(LON_profile1,LAT_profile1,DEP_profile1);
Vprofile1(isnan(Vprofile1)) = NaN; 

fig = figure;
set(gcf,'Position',[200,200,1500,700],'Color','w')

surface(LON_profile1,LAT_profile1,DEP_profile1,Vprofile1,'EdgeColor','none'); 
xlabel('Longitude (deg)');
ylabel('Latitude (deg)');
zlabel('Depth (km)')
set(gca,'FontSize',16)
set(gca,'Zdir','reverse')
zlim([0 100])
shading interp
view(55,40)
grid on;


end