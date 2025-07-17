function dBinned = doBinning(itrMat, rx, ry, x, y, dx, dy,flag)
% DOBINNING  Bin trace data (itrMat) along x, using bins centered at xgrid.
%
% Inputs:
%   itrMat  : [Nt x Ntrace], each column is a trace
%   rx      : [1 x Ntrace], horizontal location of each trace
%   xgrid   : array of bin centers [1 x nx]
%   dx      : bin width (km)
%
% Output:
%   dBinned : [Nt x nx], binned (averaged) data

    [nt, nTr] = size(itrMat);
    nx = length(x);
    ny = length(y);
    
    dBinned = zeros(nt, nx, ny);
    
    if flag
        for i = 1:nx
            for j = 1:ny
                keepx = rx>=x(i)-dx & rx<=x(i)+dx;
                keepy = ry>=y(j)-dy & ry<=y(j)+dy;
                if sum(keepx)>0 && sum(keepy)>0
        
                    rfpx = rx(logical(keepx.*keepy));
                    rfpy = ry(logical(keepx.*keepy));
        
                    rftmp = [];
                    for k = 1:length(rfpx)
                        a1 = abs(rfpx(k)-rx)<0.1;
                        a2 = abs(rfpy(k)-ry)<0.1;
        
                        nrf = find(rx==rx(logical(a1.*a2)));
        
                        rftmp = [rftmp nrf];
        
                    end
        
                    tmp = mean(itrMat(:,rftmp),2);
        
                    dBinned(:,i,j) = tmp;
                end
        
            end
        
        end
    else
        [X_grid, Y_grid] = meshgrid(x, y);
        x_coords = rx;
        y_coords = ry;
        F = scatteredInterpolant(x_coords, y_coords, zeros(size(x_coords)), 'linear', 'none');
  
        for t = 1:size(itrMat,1)

            current_data = itrMat(t,:);
            F.Values = current_data(:);
            interpolated = F(X_grid, Y_grid);
            dBinned(t, :, :) = permute(interpolated, [2, 1]);

            % 显示进度
            if mod(t,100) == 0
                fprintf('已完成 %.1f%% (%d/%d)\n', t/size(itrMat,1)*100, t, size(itrMat,1));
            end
        end


    end

end