function [rx, ry] = latlonToProjectedCoords(stlo, stla, gridStruct)
    % convert the coordinates from lat, lon, to projected Cartesian coordinates
    [stationX, stationY] = latlon2xy(stlo, stla, gridStruct.originLon, gridStruct.originLat);
    
    % 执行PCA分析
    coords = [stationX(:), stationY(:)];
    projection_on_principal_axis = coords * gridStruct.coeff(:, 1);  % 主轴方向投影
    projection_on_secondary_axis = coords * gridStruct.coeff(:, 2);  % 次轴方向投影

    % 台站投影位置
    rx = projection_on_principal_axis; % 主轴上的投影
    ry = projection_on_secondary_axis; % 次轴上的投影

end