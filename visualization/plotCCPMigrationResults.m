function plotCCPMigrationResults(ccpResults,migResults, gridStruct,smoothLength)
if nargin < 4
    smoothLength = 0;
end
x = gridStruct.x;
z = gridStruct.z;
stackedImage = zeros(size(ccpResults(1).img));
totalCount = 0;

for n=1:length(ccpResults)
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
for n=1:length(migResults)
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
for n=1:length(migResults)
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
subplot(311)
imagesc(x,z,d2dccp); hold on;
axis([xmin xmax 0 100])
xlabel('Distance (km)');
ylabel('Depth (km)');
title('CCP image')
set(gca,'fontsize',14)
cmax=3*rms(d2dccp(:));
caxis([-cmax cmax]);
colorbar
text(-0.12,0.98,'a)','Units','normalized','FontSize',18)

subplot(312)
imagesc(x,z,d2dmig); hold on;
axis([xmin xmax 0 100])
xlabel('Distance (km)');
ylabel('Depth (km)');
title('Migration image')
set(gca,'fontsize',14)
cmax=3*rms(d2dmig(:));
caxis([-cmax cmax]);
colorbar
text(-0.12,0.98,'b)','Units','normalized','FontSize',18)

subplot(313)
imagesc(x,z, d2dlsm); hold on;
axis([xmin xmax 0 100])
xlabel('Distance (km)');
ylabel('Depth (km)');
title('LSM image')
set(gca,'fontsize',14)
cmax=3*rms(d2dlsm(:));
caxis([-cmax cmax]);
colorbar
colormap(flipud(cmap.roma));
text(-0.12,0.98,'c)','Units','normalized','FontSize',18)

% export_fig(gcf,'./figures/ccp_mig_lsm.png','-r150')
