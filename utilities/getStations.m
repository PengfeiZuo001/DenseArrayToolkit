function stationList = getStations(DataStruct)
stationList = {};
for n = 1:length(DataStruct)
    stationList{end+1} = DataStruct(n).StationInfo.sta;
end
stationList = unique(stationList)';