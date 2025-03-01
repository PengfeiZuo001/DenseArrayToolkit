% 加载CCP结果数据
load ccp_QB.mat;
% 加载颜色映射表
load roma.mat;

% 是否选择剖面
select_profile = 0;
if select_profile
    % 选择剖面
    figure;
    set(gcf,'Position',[50 50 800 800],'Color','w')
    % 绘制三维切片图
    h = slice(X,Y,Z,V,90,90,40);
    xlabel('X (km)');
    ylabel('Y (km)');
    zlabel('Z (km)');
    set(h(:),'EdgeColor','none')
    set(gca,'ZDir','reverse')
    colormap(flipud(roma));
    cmax = rms(V(:));
    caxis([-cmax cmax]);
    view(0,90)
    axis equal
    % 通过鼠标点击选择剖面点
    [xpoints,ypoints] = ginput;

    % 生成剖面
    profileAll={};
    for n = 1:length(xpoints)-1
        profile=[xpoints(n),ypoints(n);xpoints(n+1),ypoints(n+1)];
        profileAll{n} = profile;
    end
end

% 绘制CCP剖面图
plotCCPXsectionCartesian(X,Y,Z,V,gridStruct,profileAll)
