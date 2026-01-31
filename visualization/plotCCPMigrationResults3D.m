function plotCCPMigrationResults3D(ccpResults,migResults, gridStruct,smoothLength,dem)
    
    load('../visualization/colormap/roma.mat')

    if ~isempty(dem)
        lonp = gridStruct.principalAxisLatLon(:,1);
        latp = gridStruct.principalAxisLatLon(:,2);
        elevp = interp2(dem.demLon,dem.demLat,dem.demZ,lonp,latp)/1000.;
    end

    x = gridStruct.x;
    y = gridStruct.y;
    z = gridStruct.z;
    [Z,X,Y] = ndgrid(z,x,y);
    Z = Z(:);   X = X(:);   Y = Y(:);
    zmax = max(z);
    xmin = 0;
    xmax = max(x);

    ccpImage = ccpResults.img;
    migImage = migResults.mig;
    lsmImage = migResults.migls;

    Fi1 = scatteredInterpolant(X,Y,Z,ccpImage(:));
    Fi2 = scatteredInterpolant(X,Y,Z,migImage(:));
    Fi3 = scatteredInterpolant(X,Y,Z,lsmImage(:));
    
    % all stations
    xx = gridStruct.rx(:,1);
    yy = gridStruct.ry(:,2);

    xprof = [];
    Vprofile_sum1 = [];
    Vprofile_sum2 = [];
    Vprofile_sum3 = [];
    npt = length(xx);
    dist = zeros(1,npt-1);
    nz = length(z);

    for i = 1:npt-1
        ninp = 1;
        xxi = linspace(xx(i),xx(i+1),ninp);
        yyi = linspace(yy(i),yy(i+1),ninp);

        x_pro = repmat(xxi,nz,1);
        y_pro = repmat(yyi,nz,1);
        depth_pro = repmat(z',1,ninp);

        Vprofile1 = Fi1(x_pro,y_pro,depth_pro);
        Vprofile1(isnan(Vprofile1)) = NaN;
        Vprofile_sum1 = [Vprofile_sum1 Vprofile1];

        Vprofile2 = Fi2(x_pro,y_pro,depth_pro);
        Vprofile2(isnan(Vprofile2)) = NaN;
        Vprofile_sum2 = [Vprofile_sum2 Vprofile2];

        Vprofile3 = Fi3(x_pro,y_pro,depth_pro);
        Vprofile3(isnan(Vprofile3)) = NaN;
        Vprofile_sum3 = [Vprofile_sum3 Vprofile3];

        dist(1,i) = sqrt((xx(i+1)-xx(i))^2 + (yy(i+1)-yy(i))^2 );
        if i >1
            xprof = [xprof linspace(0,dist(1,i),ninp)+(xprof(end))];
        elseif i==1
            xprof = linspace(0,dist(1,i),ninp);
        end

    end

    ngrid = ninp + smoothLength;


    cmap = load('../visualization/colormap/roma.mat');
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
    imagesc(xprof,z,Vprofile_sum1); hold on;
    axis([xmin xmax 0 zmax])
    % xlabel('Distance (km)');
    ylabel('Depth (km)');
    title('CCP image')
    set(gca,'fontsize',14)
    cmax=2*rms(Vprofile_sum1(:));
    caxis([-cmax cmax]);
    % colorbar
    set(gca,'XTickLabel',[])
    text(-0.12,0.98,'a)','Units','normalized','FontSize',18)
    % axis equal

    subplot(10,1,5:7)
    imagesc(xprof,z,Vprofile_sum2); hold on;
    axis([xmin xmax 0 zmax])
    % xlabel('Distance (km)');
    ylabel('Depth (km)');
    title('Migration image')
    set(gca,'fontsize',14)
    cmax=2*rms(Vprofile_sum2(:));
    caxis([-cmax cmax]);
    % colorbar
    set(gca,'XTickLabel',[])
    text(-0.12,0.98,'b)','Units','normalized','FontSize',18)
    % axis equal

    subplot(10,1,8:10)
    imagesc(xprof,z, Vprofile_sum3); hold on;
    axis([xmin xmax 0 zmax])
    xlabel('Distance (km)');
    ylabel('Depth (km)');
    title('LSM image')
    set(gca,'fontsize',14)
    cmax=2*rms(Vprofile_sum3(:));
    caxis([-cmax cmax]);
    % colorbar
    colormap(flipud(cmap.roma));
    text(-0.12,0.98,'c)','Units','normalized','FontSize',18)


end
