function plotCCPMigrationResults(dimg,dmig,dmigls,x,z)
d2dccp=mean(dimg,3);
d2dmig=mean(dmig,3);
d2dlsm=mean(dmigls,3);
xmin = min(x);
xmax = max(x);
% load colormap
cmap = load('./visualization/colormap/roma.mat');
figure();
set(gcf,'Position',[100 100 800 1200],'color','w')
subplot(311)
imagesc(x,z,d2dccp/max(abs(d2dccp(:)))); hold on;
axis([xmin xmax 0 100])
xlabel('Distance (km)');
ylabel('Depth (km)');
title('CCP image')
set(gca,'fontsize',14)
cmax=0.3;
caxis([-cmax cmax]);
colorbar
text(-0.12,0.98,'a)','Units','normalized','FontSize',18)

subplot(312)
imagesc(x,z,d2dmig/max(abs(d2dmig(:)))); hold on;
axis([xmin xmax 0 100])
xlabel('Distance (km)');
ylabel('Depth (km)');
title('Migration image')
set(gca,'fontsize',14)
caxis([-cmax cmax]);
colorbar
text(-0.12,0.98,'b)','Units','normalized','FontSize',18)

subplot(313)
imagesc(x,z, d2dlsm/max(abs(d2dlsm(:)))); hold on;
axis([xmin xmax 0 100])
xlabel('Distance (km)');
ylabel('Depth (km)');
title('LSM image')
set(gca,'fontsize',14)
caxis([-cmax cmax]);
colorbar
colormap(flipud(cmap.roma));
text(-0.12,0.98,'c)','Units','normalized','FontSize',18)

export_fig(gcf,'./figures/ccp_mig_lsm.png','-r150')
