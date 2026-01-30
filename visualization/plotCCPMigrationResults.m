function plotCCPMigrationResults(ccpResults,migResults, gridStruct,smoothLength,dem)
if nargin < 4
    smoothLength = 0;
end
if nargin < 5
    dem = [];
end

% extract surface elevation from DEM file if provided
if ~isempty(dem)
    lonp = gridStruct.principalAxisLatLon(:,1);
    latp = gridStruct.principalAxisLatLon(:,2);
    elevp = interp2(dem.demLon,dem.demLat,dem.demZ,lonp,latp)/1000.;
end

x = gridStruct.x;
z = gridStruct.z;
zmax = max(z);
stackedImage = zeros(size(ccpResults(1).img));
totalCount = 0;

for n=1:length(ccpResults)-1
    img = ccpResults(n).img;
    count = ccpResults(n).count;
    if smoothLength>0
        K = (1/smoothLength^2)*ones(smoothLength,smoothLength);
        img = conv2(img,K,'same');
        count = conv2(count,K,'same');
    end
    stackedImage = stackedImage+img;
    totalCount = totalCount+count;
end
d2dccp = stackedImage./max(totalCount,1);
% d2dccp = V/length(ccpResults);

stackedImage = zeros(size(migResults(1).mig));
for n=1:length(migResults)-1
    mig=migResults(n).mig;
    if smoothLength>0
        K = (1/smoothLength^2)*ones(smoothLength,smoothLength);
        mig = conv2(mig,K,'same');
    end
    stackedImage = stackedImage+mig;
end

d2dmig = stackedImage/length(migResults);
% d2dmig = V-mean(V(:));

stackedImage = zeros(size(migResults(1).migls));
for n=1:length(migResults)-1
    migls=migResults(n).migls;
    if smoothLength>0
        K = (1/smoothLength^2)*ones(smoothLength,smoothLength);
        migls = conv2(migls,K,'same');
    end
    stackedImage = stackedImage+migls;
end

d2dlsm = stackedImage/length(migResults);
% d2dlsm = V-mean(V(:));

xmin = 0;
xmax = max(x)+min(x);

% load colormap
cmap = load('./visualization/colormap/roma.mat');
figure();
set(gcf,'Position',[100 100 800 1200],'color','w')
if ~isempty(dem)
    subplot(10,1,1)
    plot(x,elevp,'k','linewidth',1); hold on;
    patch([x(1);x(:);x(end)],[0;elevp;0],[0.8 0.8 0.8]);
    ylabel('Elev. (km)')
    xlim([0 xmax])
    ylim([2 5])
    set(gca,'XTickLabel',[])
    set(gca,'fontsize',12)
end
subplot(10,1,2:4)
imagesc(x,z,d2dccp); hold on;
axis([xmin xmax 0 zmax])
% xlabel('Distance (km)');
ylabel('Depth (km)');
title('CCP image')
set(gca,'fontsize',14)
cmax=2*rms(d2dccp(:));
caxis([-cmax cmax]);
% colorbar
set(gca,'XTickLabel',[])
text(-0.12,0.98,'a)','Units','normalized','FontSize',18)
% axis equal

subplot(10,1,5:7)
imagesc(x,z,d2dmig); hold on;
axis([xmin xmax 0 zmax])
% xlabel('Distance (km)');
ylabel('Depth (km)');
title('Migration image')
set(gca,'fontsize',14)
cmax=2*rms(d2dmig(:));
caxis([-cmax cmax]);
% colorbar
set(gca,'XTickLabel',[])
text(-0.12,0.98,'b)','Units','normalized','FontSize',18)
% axis equal

subplot(10,1,8:10)
imagesc(x,z, d2dlsm); hold on;
axis([xmin xmax 0 zmax])
xlabel('Distance (km)');
ylabel('Depth (km)');
title('LSM image')
set(gca,'fontsize',14)
cmax=2*rms(d2dlsm(:));
caxis([-cmax cmax]);
% colorbar
colormap(flipud(cmap.roma));
text(-0.12,0.98,'c)','Units','normalized','FontSize',18)
% axis equal

% export_fig(gcf,'./figures/ccp_mig_lsm.png','-r150')
