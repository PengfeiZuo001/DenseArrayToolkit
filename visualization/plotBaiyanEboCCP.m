load BaiyanEbo_ccp.mat
load roma.mat;
% Baiyan Ebo location
lon = 109.9250;
lat = 41.7666;

[x,y]=latlon2xy(lon,lat,gridStruct.originLon,gridStruct.originLat);
figure;
h = slice(X,Y,Z,V_smooth,x,y,40);
xlabel('X (km)');
ylabel('Y (km)');
zlabel('Z (km)');
set(h(:),'EdgeColor','none')
set(gca,'ZDir','reverse')
colormap(flipud(roma));
cmax = rms(V_smooth(:));
caxis([-cmax cmax]);