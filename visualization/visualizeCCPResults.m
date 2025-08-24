function [profileStruct] = visualizeCCPResults(ccpResult, gridStruct, options)
%% visualizeCCPResults - Visualize CCP stacking results with structure-oriented filtering
% This function creates visualizations of CCP stacking results, including 3D volume
% visualization and 2D profiles with optional structure-oriented filtering. It uses
% plotCCPXsectionCartesian for profile display and adds structure-oriented filtering
% capabilities.
%
% Usage:
%   profileStruct = visualizeCCPResults(ccpResult, gridStruct)
%   profileStruct = visualizeCCPResults(ccpResult, gridStruct, options)
%
% Inputs:
%   ccpResult  - Structure containing CCP stacking results
%                .X, .Y, .Z: Coordinate grids
%                .V: Image volume (3D array)
%   gridStruct - Structure containing velocity model and coordinate information
%   options    - (Optional) Structure containing visualization parameters
%                .profileType: 'interactive' (default) or 'predefined'
%                             In interactive mode, left-click to add points,
%                             right-click to finish selection
%                .smoothingParams: Structure for filtering parameters
%                    .radius: Smoothing radius (default: 3)
%                    .eps: Regularization parameter (default: 0.01)
%                    .order: Smoothing order (default: 2)
%
% Outputs:
%   profileStruct - Structure containing profile information:
%                   .distance: Distance along profile
%                   .depth: Depth points
%                   .amplitude: Original profile data
%                   .filtered: Structure-oriented filtered data
%                   .dip: Calculated structural dips
%
% The function provides interactive profile selection with real-time visual
% feedback and integrates with plotCCPXsectionCartesian for consistent
% profile visualization across the toolkit.

%% Input parameter checking and defaults
if nargin < 3
    options = struct();
end

if ~isfield(options, 'profileType')
    options.profileType = 'interactive';
end

if ~isfield(options, 'displayMode')
    options.displayMode = 'both';
end

% if ~isfield(options, 'smoothingParams')
%     options.smoothingParams = struct('radius', 3, 'eps', 0.01, 'order', 2);
% end

%% Load required data
% Load colormap
cmap = load('./visualization/colormap/roma.mat');
%% Prepare visualization data
% Extract coordinates
X = ccpResult.X;
Y = ccpResult.Y;
Z = ccpResult.Z;
V = ccpResult.V;

%% 3D Volume Visualization

figure('Position', [0 0 1000 1000], 'Color', 'w');
ax1 = subplot(1,1,1);
hold on;

% Plot volume slices
h = slice(X, Y, Z, V, [], [], 100);
colormap(ax1, flipud(cmap.roma));
cmax = 2*rms(V(:));
caxis([-cmax, cmax]);

% Set axes properties
xlabel('X (km)');
ylabel('Y (km)');
zlabel('Z (km)');
set(h(:), 'EdgeColor', 'none');
set(ax1, 'ZDir', 'reverse', 'YDir', 'reverse', 'FontSize', 16);

% Plot station locations
rx = gridStruct.rx(:,1);
ry = gridStruct.ry(:,2);
scatter(ax1, rx, ry, 100, '^', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'k');

% Set DEM display properties
zlim([ax1.ZLim]);
xlim([ax1.XLim]);
ylim([ax1.YLim]);
axis equal
axis xy


%% Profile Selection and Processing
depth0 = 0:0.5:max(Z(:));
if strcmp(options.profileType, 'interactive')
    figure(gcf);
    view(0,90);
    hold on;
    
    % Initialize arrays for points and lines
    xpoints = [];
    ypoints = [];
    h_points = [];  % Handle array for point markers
    h_lines = [];   % Handle array for connecting lines
    
    % Interactive point selection with visual feedback
    while true
        [x, y, button] = ginput(1);
        
        % Right click ends selection
        if button > 1
            break;
        end
        
        % Store point coordinates
        xpoints = [xpoints; x];
        ypoints = [ypoints; y];
        
        % Plot new point
        h = plot(ax1, x, y, 'ro', 'MarkerSize', 10, 'LineWidth', 2);
        h_points = [h_points; h];
        
        % Draw line connecting points
        if length(xpoints) > 1
            h_line = plot(ax1, [xpoints(end-1), x], [ypoints(end-1), y], 'r-', 'LineWidth', 2);
            h_lines = [h_lines; h_line];
        end
    end
    
    % Ensure at least 2 points were selected
    if length(xpoints) < 2
        error('At least 2 points must be selected to define a profile');
    end
    
    hold off;
else
    % Use predefined profile points if provided
    if isfield(options, 'profilePoints')
        xpoints = options.profilePoints(:,1);
        ypoints = options.profilePoints(:,2);
        h = plot(ax1, xpoints, ypoints, 'ro', 'MarkerSize', 10, 'LineWidth', 2);
        h_line = plot(ax1, xpoints, ypoints, 'r-', 'LineWidth', 2);
    else
        error('Profile points must be provided for predefined profile type');
    end
end

%% Process Profiles
% Convert selected points to profile format for plotCCPXsectionCartesian
profileAll = {};
for n = 1:length(xpoints)-1
    profile = [xpoints(n), ypoints(n); xpoints(n+1), ypoints(n+1)];
    profileAll{n} = profile;
end

% Call plotCCPXsectionCartesian for profile visualization
[distAll, depthAll, VAll] = plotCCPXsectionCartesian(X, Y, Z, V, gridStruct, profileAll);

% Calculate structural dips and apply filtering if requested
if isfield(options, 'smoothingParams')
    % Normalize data for dip calculation
    VAll_norm = VAll/max(abs(VAll(:)));
    
    % Calculate structural dips
    dip = str_dip2d(VAll_norm);
    
    % Apply structure-oriented smoothing
    VAll_filtered = str_pwsmooth_lop2d(VAll, dip, ...
                                     options.smoothingParams.radius, ...
                                     options.smoothingParams.order, ...
                                     options.smoothingParams.eps);
                                     
    % Display filtered results
    figure;
    set(gcf,'Position',[0 0 1000 1000], 'Color', 'w')
    subplot(211)
    imagesc(distAll(1,:), depthAll(:,1), VAll);
    title('Original Profile');
    colormap(flipud(cmap.roma));
    cmax = 2*rms(VAll(:));
    caxis([-cmax, cmax]);
    colorbar;
    xlabel('Distance (km)');
    ylabel('Depth (km)');
    set(gca, 'FontSize', 16);

    subplot(212)
    imagesc(distAll(1,:), depthAll(:,1), VAll_filtered);
    title('Structure-Oriented Filtered Profile');
    colormap(flipud(cmap.roma));
    caxis([-cmax, cmax]);
    colorbar;
    xlabel('Distance (km)');
    ylabel('Depth (km)');
    set(gca, 'FontSize', 16);

    % Store results in profile structure
    profileStruct = struct('distance', distAll, ...
        'depth', depthAll, ...
        'amplitude', VAll, ...
        'filtered', VAll_filtered, ...
        'dip', dip);
else

    % Store results in profile structure
    profileStruct = struct('distance', distAll, ...
        'depth', depthAll);
end
