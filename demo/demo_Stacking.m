%% Main Function for Receiver Function Computation and Visualization
% This script demonstrates a typical workflow from reading data to computing receiver functions, 
% and then visualizing the stacked results.
% The workflow includes the following steps:
%   0. Setup paths and parameters
%   1. Read data
%   2. Preprocessing
%   3. Compute receiver functions
%   4. Stack receiver functions
%   5. Visualization

clear; clc; close all;

%% 0. Setup Paths and Parameters
% 0.1 Load configuration
%   Uses the custom function loadConfig() to load configuration files, 
%   such as data paths and processing parameters.
config = loadConfig();

% 0.2 Use configuration parameters
%   Extracts important parameters from the configuration file as local variables
%   for use in subsequent functions.
dataFolder         = config.dataFolder;
PreprocessingParam = config.PreprocessingParam;
MigParam           = config.MigParam;        % Can be ignored or removed if not needed in this script
DeconvParam        = config.DeconvParam;

%% 1. Read Data
%   Uses read_SAC() function to read SAC-format seismic data from the dataFolder directory,
%   and encapsulates the data into the DataStruct structure.
DataStruct = read_SAC(dataFolder);

%% 2. Preprocessing
%   Calls preprocessing() function to perform preprocessing operations on the data:
%   For example, remove mean, detrend, bandpass filtering, time window cutting, etc.
%   (Specific operations are determined by PreprocessingParam).
DataStruct = preprocessing(DataStruct, PreprocessingParam);

%% 3. Compute Receiver Functions
%   Calls deconv() function to compute receiver functions from seismic records 
%   using inverse filtering methods.
%   Different deconvolution methods, water level factors, and smoothing parameters
%   can be configured in DeconvParam.
DataStruct = deconv(DataStruct, DeconvParam);

%% 4. Stack Receiver Functions
%   Stacks receiver functions from all events for the same station to obtain
%   more stable and clearer results.
%   stackCommonStationGather() returns the stacked data seisout, depth coordinates depth0,
%   and optional Moho information mohoStruct.
[seisout, depth0, mohoStruct] = stackCommonStationGather(DataStruct);

%% 5. Visualization
%   After computing and stacking receiver functions, visualize and inspect the results.

% 5.1 Visualize example single seismic trace
trace_index = 100;  % Specify the index of the trace to plot
%   plotWaveforms() plots the waveform of a specific seismic record in time or frequency domain
plotWaveforms(DataStruct, trace_index);

% 5.2 Visualize receiver functions for a specific station
station = getStations(DataStruct);
stationList = {station.sta};
station = stationList{1};  % Select the first station for demonstration
plotCommonStationGather(DataStruct, station);

% 5.3 Visualize receiver functions for all events
%   First get the event list, then loop through each event to plot and save as PNG files.
event     = getEvents(DataStruct);
eventList = {event.evid};
for iEvent = 1:length(eventList)
    EventID = eventList{iEvent};
    
    % Extract seismic records corresponding to this event
    gather = getCommonEventGather(DataStruct, EventID);
    
    % Use plotCommonEventGather() for plotting
    %   If 'trace' is specified, it plots each trace waveform; 
    %   other options include 'wiggle', 'stack', etc.
    plotCommonEventGather(DataStruct, EventID, 'trace','wigb');
    
    % Save current figure as PNG file with 150 dpi resolution
    % export_fig(['./figures/', EventID, '.png'], '-r150');
    
    % Close all current figures to prevent generating too many windows during batch processing
    close all
end

% 5.4 Visualize station and event distribution on map
%   plotStations() plots station locations based on station information in DataStruct,
%   and passes external map file to function for background plotting.
plotStations(DataStruct, 'Baiyanebo_DEM.mat');

%   plotEvents() plots event location distribution
plotEvents(DataStruct);
