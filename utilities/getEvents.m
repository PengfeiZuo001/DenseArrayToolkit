function eventList = getEvents(DataStruct)
eventList = {};
for n = 1:length(DataStruct)
    eventList{end+1} = DataStruct(n).EventInfo.evid;
end
eventList = sort(unique(eventList)');