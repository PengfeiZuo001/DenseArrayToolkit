function [stackedResult] = stackImagingResults(resultsArray)
%% stackImagingResults - Stack multiple imaging results into a single volume
% This function combines multiple migration or CCP imaging results by 
% averaging or summing them to create a single stacked image volume.
% For migration results, it averages the least-squares migration images.
% For CCP results, it sums the images and normalizes by the count.
%
% Usage:
%   stackedResult = stackImagingResults(resultsArray)
%
% Inputs:
%   resultsArray - Array of imaging result structures containing either:
%                  - Migration results: each element has .migls field (3D array)
%                  - CCP results: each element has .img and .count fields
%
% Outputs:
%   stackedResult - Structure containing stacked imaging results:
%                   .V: Stacked image volume (3D array)
%                   .X, .Y, .Z: Coordinate grids (copied from first result)
%                   .count: Total count for CCP stacking (CCP results only)
%
% Example:
%   % Stack migration results
%   migImage = stackImagingResults(migResults);
%
%   % Stack CCP results  
%   ccpImage = stackImagingResults(ccpResults);

%% Input validation
if isempty(resultsArray)
    error('resultsArray cannot be empty');
end

if ~isstruct(resultsArray)
    error('resultsArray must be a structure array');
end

%% Determine result type and stack accordingly
if isfield(resultsArray(1), 'migls')
    % Migration results - average all migls fields
    stackedResult = stackMigrationResults(resultsArray);
elseif isfield(resultsArray(1), 'img') && isfield(resultsArray(1), 'count')
    % CCP results - sum images and normalize by count
    stackedResult = stackCCPResults(resultsArray);
else
    error('Unknown result type. Expected migration (.migls) or CCP (.img, .count) fields');
end

end

%% Helper function for stacking migration results
function [stacked] = stackMigrationResults(resultsArray)
% Stack migration results by averaging least-squares migration images
    
    % Initialize with first result
    firstResult = resultsArray(1);
    stackedmig = zeros(size(firstResult.mig));
    stackedmigls = zeros(size(firstResult.migls));
    
    % Sum all migration images
    for i = 1:length(resultsArray)
        stackedmig = stackedmigls + resultsArray(i).mig;
        stackedmigls = stackedmigls + resultsArray(i).migls;
    end
    
    % Average the stacked volume
    stackedmigls = stackedmigls ./ length(resultsArray);
    
    % Create output structure
    stacked = struct();
    stacked.mig = stackedmig;
    stacked.migls = stackedmigls;
    
    % Copy coordinate grids if available
    if isfield(firstResult, 'x')
        stacked.x = firstResult.x;
    end
    if isfield(firstResult, 'y')
        stacked.y = firstResult.y;
    end
    if isfield(firstResult, 'z')
        stacked.z = firstResult.z;
    end
end

%% Helper function for stacking CCP results
function [stacked] = stackCCPResults(resultsArray)
% Stack CCP results by summing images and normalizing by count
    
    % Initialize with first result
    firstResult = resultsArray(1);
    stackedVolume = zeros(size(firstResult.img));
    totalCount = zeros(size(stackedVolume));
    
    % Sum all images and counts
    for i = 1:length(resultsArray)
        stackedVolume = stackedVolume + resultsArray(i).img;
        totalCount = totalCount + resultsArray(i).count;
    end

    % Normalize by total count (avoid division by zero)
    Vn = zeros(size(stackedVolume));
    mask = totalCount > 0;
    Vn(mask) = stackedVolume(mask) ./ totalCount(mask);
    Vn(~mask) = 0;   % or Nan

    % Create output structure
    stacked = struct();
    stacked.img = Vn;
    stacked.count = totalCount;
    
    % Copy coordinate grids if available
    if isfield(firstResult, 'X')
        stacked.X = firstResult.X;
    end
    if isfield(firstResult, 'Y')
        stacked.Y = firstResult.Y;
    end
    if isfield(firstResult, 'Z')
        stacked.Z = firstResult.Z;
    end
end
