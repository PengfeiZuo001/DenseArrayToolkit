# 定义区域
region=89/98/36/42
grd_file=earth_relief_15s.grd
# 下载高分辨率的地表高程数据
gmt grdcut @earth_relief_15s -R$region -G$grd_file
