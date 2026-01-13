function [out] = radon3d_op(in, Param, operator)
% RADON3D_OP - 3D Radon transform for irregular station layouts
%
% This function extends the standard 3D Radon transform to handle non-regular
% grid configurations, which is essential for real-world seismic arrays where
% stations are not arranged in a perfect rectangular grid.
%
% INPUTS:
%   in        - Input data array (model space or data space)
%               For operator = 1: model [nt, npx, npy] or [nt, nv]
%               For operator = -1: data [nt, nhx] (non-regular) or [nt, nhx, nhy] (regular)
%   Param     - Parameter structure containing:
%               .hx    : x-offset vector [nhx x 1] (receiver x-coordinates)
%               .hy    : y-offset vector [nhy x 1] (receiver y-coordinates)
%               .nt    : Number of time samples
%               .dt    : Time sampling interval (seconds)
%               .type  : Moveout type (1=linear, 2=parabolic, 3=hyperbolic)
%               .px    : x-slowness parameters [npx x 1] (types 1-2)
%               .py    : y-slowness parameters [npy x 1] (types 1-2)
%               .v     : Velocity parameters [nv x 1] (type 3 only)
%               .isGrid: Grid type flag (1 = regular grid, 0 = non-regular stations)
%   operator  - Operation direction:
%               1  = Forward operator (model -> data)
%               -1 = Adjoint operator (data -> model)
%
% OUTPUT:
%   out       - Output data array (data space or model space)
%               For operator = 1 and isGrid=1: data [nt, nhx, nhy]
%               For operator = 1 and isGrid=0: data [nt, nhx]
%               For operator = -1: model [nt, npx, npy] or [nt, nv]
%
% KEY FEATURES:
%   - Supports both regular grid and non-regular station layouts
%   - For non-regular grids, stations are treated as 1D array with combined (hx, hy) pairs
%   - Maintains adjoint property for iterative inversion algorithms
%
% Author: MATLAB DenseArrayToolkit
% References: Adapted from regular grid Radon transform for seismic applications

% Extract parameters from structure
hx = Param.hx;      % x-offset coordinates of receivers
hy = Param.hy;      % y-offset coordinates of receivers
nt = Param.nt;      % Number of time samples
dt = Param.dt;      % Time sampling interval (seconds)
type = Param.type;  % Transform type (1=linear, 2=parabolic, 3=hyperbolic)
isGrid = Param.isGrid;  % Grid type flag (1=regular, 0=non-regular)

% Calculate array dimensions
nhx = length(hx);   % Number of x-offsets
nhy = length(hy);   % Number of y-offsets

% =========================================================================
% INITIALIZATION: HANDLE BOTH REGULAR AND NON-REGULAR GRIDS
% =========================================================================

if operator == 1      % FORWARD OPERATOR: Model -> Data
    if type == 3      % Hyperbolic transform (velocity domain)
        nv = length(Param.v);           % Number of velocity parameters
        if isGrid
            d = zeros(nt, nhx, nhy);    % Regular grid: 3D data array
        else
            d = zeros(nt, nhx);         % Non-regular: 2D data array (stations only)
        end
        m = in;                         % Input model: [nt, nv]
    else              % Linear or Parabolic transform (slowness domain)
        npx = length(Param.px);         % Number of x-slowness parameters
        npy = length(Param.py);         % Number of y-slowness parameters
        if isGrid
            d = zeros(nt, nhx, nhy);    % Regular grid: 3D data array
        else
            d = zeros(nt, nhx);         % Non-regular: 2D data array
        end
        m = in;                         % Input model: [nt, npx, npy]
    end
    
else                % ADJOINT OPERATOR: Data -> Model
    d = in;  % Input data (shape depends on isGrid flag)
    
    if type == 3    % Hyperbolic transform (velocity domain)
        nv = length(Param.v);           % Number of velocity parameters
        m = zeros(nt, nv);              % Output model: [nt, nv]
    else            % Linear or Parabolic transform (slowness domain)
        npx = length(Param.px);         % Number of x-slowness parameters
        npy = length(Param.py);         % Number of y-slowness parameters
        m = zeros(nt, npx, npy);        % Output model: [nt, npx, npy]
    end
end

% =========================================================================
% RADON TRANSFORM IMPLEMENTATION WITH GRID TYPE SUPPORT
% =========================================================================

switch type
    case 1  % LINEAR RADON TRANSFORM
        px = Param.px;  % x-slowness parameters
        py = Param.py;  % y-slowness parameters
        
        for itau = 1:nt           % Loop over intercept time
            for ipx = 1:length(px) % Loop over x-slowness parameters
                for ipy = 1:length(py) % Loop over y-slowness parameters
                    
                    if isGrid
                        % REGULAR GRID: Full 2D spatial sampling
                        for ihx = 1:nhx    % Loop over x-offsets
                            for ihy = 1:nhy % Loop over y-offsets
                                % Compute travel time for regular grid
                                t = (itau-1)*dt + px(ipx)*hx(ihx) + py(ipy)*hy(ihy);
                                it = floor(t/dt) + 1;  % Convert to sample index
                                
                                if it >= 1 && it <= nt
                                    if operator == 1
                                        % Forward: Map model to 3D data grid
                                        d(it, ihx, ihy) = d(it, ihx, ihy) + m(itau, ipx, ipy);
                                    else
                                        % Adjoint: Map 3D data grid to model
                                        m(itau, ipx, ipy) = m(itau, ipx, ipy) + d(it, ihx, ihy);
                                    end
                                end
                            end
                        end
                    else
                        % NON-REGULAR GRID: Stations treated as 1D array
                        % Each station has unique (hx, hy) coordinates
                        for ihx = 1:nhx    % Loop over stations
                            % For non-regular grid, hx and hy arrays have same length
                            % Each index corresponds to one station with (hx, hy) coordinates
                            t = (itau-1)*dt + px(ipx)*hx(ihx) + py(ipy)*hy(ihx);
                            it = floor(t/dt) + 1;  % Convert to sample index
                            
                            if it >= 1 && it <= nt
                                if operator == 1
                                    % Forward: Map model to 2D station data
                                    d(it, ihx) = d(it, ihx) + m(itau, ipx, ipy);
                                else
                                    % Adjoint: Map 2D station data to model
                                    m(itau, ipx, ipy) = m(itau, ipx, ipy) + d(it, ihx);
                                end
                            end
                        end
                    end
                end
            end
        end
        
    case 2  % PARABOLIC RADON TRANSFORM
        % Note: Current implementation only supports regular grid for parabolic
        % This is a limitation that should be addressed in future versions
        px = Param.px;  % x-curvature parameters
        py = Param.py;  % y-curvature parameters
        
        for itau = 1:nt           % Loop over intercept time
            for ipx = 1:length(px) % Loop over x-curvature parameters
                for ipy = 1:length(py) % Loop over y-curvature parameters
                    for ihx = 1:nhx    % Loop over x-offsets
                        for ihy = 1:nhy % Loop over y-offsets
                            % Compute travel time with quadratic dependence
                            t = (itau-1)*dt + px(ipx)*hx(ihx)^2 + py(ipy)*hy(ihy)^2;
                            it = floor(t/dt) + 1;  % Convert to sample index
                            
                            if it >= 1 && it <= nt
                                if operator == 1
                                    % Forward: Accumulate model contribution
                                    d(it, ihx, ihy) = d(it, ihx, ihy) + m(itau, ipx, ipy);
                                else
                                    % Adjoint: Accumulate data contribution
                                    m(itau, ipx, ipy) = m(itau, ipx, ipy) + d(it, ihx, ihy);
                                end
                            end
                        end
                    end
                end
            end
        end

    case 3  % HYPERBOLIC RADON TRANSFORM
        % Note: Current implementation only supports regular grid for hyperbolic
        % This is a limitation that should be addressed in future versions
        v = Param.v;  % Velocity parameters
        
        for itau = 1:nt           % Loop over zero-offset time
            for iv = 1:length(v)   % Loop over velocity parameters
                for ihx = 1:nhx    % Loop over x-offsets
                    for ihy = 1:nhy % Loop over y-offsets
                        % Compute squared offset distance
                        r2 = hx(ihx)^2 + hy(ihy)^2;
                        
                        % Compute hyperbolic travel time
                        t = sqrt(((itau-1)*dt)^2 + r2/v(iv)^2);
                        it = floor(t/dt) + 1;  % Convert to sample index
                        
                        if it >= 1 && it <= nt
                            if operator == 1
                                % Forward: Map velocity model to data
                                d(it, ihx, ihy) = d(it, ihx, ihy) + m(itau, iv);
                            else
                                % Adjoint: Map data to velocity model
                                m(itau, iv) = m(itau, iv) + d(it, ihx, ihy);
                            end
                        end
                    end
                end
            end
        end
end

% =========================================================================
% OUTPUT SELECTION
% =========================================================================
if operator == 1
    out = d;  % Forward: Return data space
else
    out = m;  % Adjoint: Return model space
end
end
