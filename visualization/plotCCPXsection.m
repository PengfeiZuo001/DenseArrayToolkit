function plotCCPXsection(ccpResult,stlo,stla,profile)
%% create 2D cross-section
centers = [cell2mat(ccpResult.CCPGrid(:,2)) cell2mat(ccpResult.CCPGrid(:,1))];
radii = km2deg(cell2mat(ccpResult.CCPGrid(:,3)));
F = ccpResult.rf;
depth0 = 0:0.5:200;
% load colormap
cmap = load('./visualization/colormap/roma.mat');
figure
set(gcf,'Position',[0 0 1600 800],'Color','w')
scatter(stlo,stla,200,'^','filled'); hold on;
for n = 1:length(profile)
    lon1 = profile{n}(1,1);
    lat1 = profile{n}(1,2);
    lon2 = profile{n}(2,1);
    lat2 = profile{n}(2,2);
    nlatlon = 100;
    [latp,lonp] = gcwaypts(lat1,lon1,lat2,lon2,nlatlon);
    [deg0,az0]= distance(lat1,lon1,latp,lonp);
    % degree to distance
    dist0 = deg0*2*pi*6371/360;


    % select the ccp bin
    keepm = [];
    for m = 1:length(ccpResult.CCPGrid)
        lons = lonp * pi/180;
        lats = latp * pi/180;
        tlon = ccpResult.CCPGrid{m,2} * pi/180;
        tlat = ccpResult.CCPGrid{m,1} * pi/180;
        dlon = lons - tlon;
        dlat = lats - tlat;
        a = (sin(dlat/2)).^2 + cos(lats) .* cos(tlat) .* (sin(dlon/2)).^2;
        angles = 2 .* atan2(sqrt(a),sqrt(1-a));
        dist = 6371 * angles;
        Indx = find(dist <= ccpResult.CCPGrid{m,3});
        if ~isempty(Indx)
            keepm = [keepm m];
        end
    end
   
    % plot the geometry of RFs
    LAT_profile = repmat(latp',length(depth0),1);
    LON_profile = repmat(lonp',length(depth0),1);
    DEP_profile = repmat(depth0',1,length(latp));
    DIST_profile = repmat(dist0',length(depth0),1);
    Vprofile = F(LON_profile,LAT_profile,DEP_profile);
    Vprofile(isnan(Vprofile)) = NaN;

    scatter(lonp,latp,50,'ko','filled');
    xlabel(['Longitude' ,char(176)])
    ylabel(['Latitude' ,char(176)])
    zlabel('Depth (km)')
    set(gca,'FontSize',16)

    viscircles(centers(keepm,:),radii(keepm)); hold on;
    % viscircles(centers,radii); hold on;
    surface(LON_profile,LAT_profile,-DEP_profile,Vprofile,'EdgeColor','none')
    zlim([-100 0])
    colormap(flipud(cmap.roma))
    caxis([-0.1 0.1])
    colorbar
    grid on; box on;
    view(140,20)
    % saveas(gcf, 'example_ccp.png','png')
end