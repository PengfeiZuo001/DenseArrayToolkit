# 定义区域
#region=89/98/36/42
#region=89/102/36/42 # Qilian area
region=102/106/28/32
grd_file=earth_relief_30s.grd
# 下载高分辨率的地表高程数据
gmt grdcut @earth_relief_30s -R$region -G$grd_file
