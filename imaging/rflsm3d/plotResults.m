function plotResults(mig,pre_rfm,lsmig,pre_rflsm,param)

    param = param.paramMig;
% paramters
    x = param.x ;
    y = param.y ;
    z = param.z;

    rx = param.rx(:,1);
    ry = param.ry(:,2);

%% plot migration images
    nnx = find(y<105 & y>95,1);    % y axis  km
    nny = find(x<155 & x>145,1);    % x axis  km

    [idx1,~] = find(rx>145 & rx<155,15);
    [idx2,~] = find(ry>98 & ry<105,15);
    
    figure
    set(gcf,'Position',[1,200,1950,700],'Color','white')
    subplot(2,5,[1 2 3])
    imagesc(x,z,mig(:,:,nnx)./max(mig(:)),'CDataMapping','scaled','Interpolation','bilinear');
    title('Migration image')
    clim([-1 1]);
    axis equal
    % colormap(flipud(roma))
    ylim([0 100])
    xlim([min(x) max(x)])
    clim([-0.5 0.5])
    colorbar
    set(gca,'XMinorTick','on','YMinorTick','on');
    set(gca,'FontSize',16)
    xlabel('X (km)')
    ylabel('Depth (km)')
    hold on
    scatter(rx(idx2),zeros(1,numel(idx2)),30,'v','filled',...
                  'MarkerFaceColor',[1 0 0]);
    
    hold off
    text(-15,-12,'(a)','FontSize',20)
    
    
    subplot(2,5,[6 7 8])
    imagesc(x,z,lsmig(:,:,nnx)./max(lsmig(:)),'CDataMapping','scaled','Interpolation','bilinear');
    title('LSM image')
    axis equal
    % colormap(flipud(roma))
    ylim([0 100])
    xlim([min(x) max(x)])
    clim([-0.5 0.5])
    colorbar
    set(gca,'XMinorTick','on','YMinorTick','on');
    set(gca,'FontSize',16)
    xlabel('X (km)')
    ylabel('Depth (km)')
    hold on
    scatter(rx(idx2),zeros(1,numel(idx2)),30,'v','filled',...
                  'MarkerFaceColor',[1 0 0]);

    hold off
    text(-15,-12,'(c)','FontSize',20)
    
    subplot(2,5,[4 5])
    imagesc(y,z,squeeze(mig(:,nny,:))./max(mig(:)),'CDataMapping','scaled','Interpolation','bilinear');
    title('Migration image')
    clim([-1 1]);
    axis equal
    % colormap(flipud(roma))
    ylim([0 100])
    xlim([min(y) max(y)])
    % xlim([0 300])
    clim([-0.5 0.5])
    colorbar
    set(gca,'XMinorTick','on','YMinorTick','on');
    set(gca,'FontSize',16)
    xlabel('Y (km)')
    ylabel('Depth (km)')
    hold on
    scatter(ry(idx1),zeros(1,numel(idx1)),30,'v','filled',...
                  'MarkerFaceColor',[1 0 0]);
    hold off
    text(-15,-12,'(b)','FontSize',20)
    
    subplot(2,5,[9 10])
    imagesc(y,z,squeeze(lsmig(:,nny,:))./max(lsmig(:)),'CDataMapping','scaled','Interpolation','bilinear');
    title('LSM image')
    axis equal
    % colormap(flipud(roma))
    ylim([0 100])
    xlim([min(y) max(y)])
    clim([-0.5 0.5])
    colorbar
    set(gca,'XMinorTick','on','YMinorTick','on');
    set(gca,'FontSize',16)
    xlabel('Y (km)')
    ylabel('Depth (km)')
    hold on
    scatter(ry(idx1),zeros(1,numel(idx1)),30,'v','filled',...
                  'MarkerFaceColor',[1 0 0]);
    hold off
    text(-15,-12,'(d)','FontSize',20)


end