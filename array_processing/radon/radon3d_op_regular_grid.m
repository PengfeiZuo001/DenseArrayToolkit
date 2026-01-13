function [out] = radon3d_op_regular_grid(in, Param, operator)
% RADON3D_OP_REGULAR_GRID - 3D Radon transform operator for seismic data processing
%
% This function implements the forward and adjoint 3D Radon transform
% operators for time-domain seismic data. It supports three types of
% moveout: linear, parabolic, and hyperbolic.
%
% INPUTS:
%   in        - Input data array (model space or data space)
%               For operator = 1: model [nt, npx, npy] or [nt, nv]
%               For operator = -1: data [nt, nhx, nhy]
%   Param     - Parameter structure containing:
%               .hx    : x-offset vector [nhx x 1] (receiver x-coordinates)
%               .hy    : y-offset vector [nhy x 1] (receiver y-coordinates)
%               .nt    : Number of time samples
%               .dt    : Time sampling interval (seconds)
%               .type  : Moveout type (1=linear, 2=parabolic, 3=hyperbolic)
%               .px    : x-slowness parameters [npx x 1] (types 1-2)
%               .py    : y-slowness parameters [npy x 1] (types 1-2)
%               .v     : Velocity parameters [nv x 1] (type 3 only)
%   operator  - Operation direction:
%               1  = Forward operator (model -> data)
%               -1 = Adjoint operator (data -> model)
%
% OUTPUT:
%   out       - Output data array (data space or model space)
%               For operator = 1: data [nt, nhx, nhy]
%               For operator = -1: model [nt, npx, npy] or [nt, nv]
%
% TRANSFORM TYPES:
%   Type 1 (Linear):   t = tau + px*hx + py*hy
%   Type 2 (Parabolic): t = tau + px*hx^2 + py*hy^2  
%   Type 3 (Hyperbolic): t = sqrt(tau^2 + (hx^2+hy^2)/v^2)
%
% Author: MATLAB DenseArrayToolkit
% References: Seismic data processing literature on Radon transforms

% Extract parameters from structure for cleaner code
hx = Param.hx;    % x-offset coordinates of receivers
hy = Param.hy;    % y-offset coordinates of receivers
nt = Param.nt;    % Number of time samples
dt = Param.dt;    % Time sampling interval (seconds)
type = Param.type;% Transform type (1=linear, 2=parabolic, 3=hyperbolic)

% Calculate array dimensions
nhx = length(hx); % Number of x-offsets
nhy = length(hy); % Number of y-offsets

% =========================================================================
% INITIALIZATION: SET UP INPUT/OUTPUT ARRAYS BASED ON OPERATOR DIRECTION
% =========================================================================

if operator == 1      % FORWARD OPERATOR: Model -> Data
    % Forward operator maps from model space to data space
    if type == 3      % Hyperbolic transform (velocity domain)
        nv = length(Param.v);           % Number of velocity parameters
        d = zeros(nt, nhx, nhy);        % Initialize output data array
        m = in;                         % Input model: [nt, nv]
    else              % Linear or Parabolic transform (slowness domain)
        npx = length(Param.px);         % Number of x-slowness parameters
        npy = length(Param.py);         % Number of y-slowness parameters
        d = zeros(nt, nhx, nhy);        % Initialize output data array
        m = in;                         % Input model: [nt, npx, npy]
    end
    
else                % ADJOINT OPERATOR: Data -> Model
    % Adjoint operator maps from data space to model space
    d = in;  % Input data: [nt, nhx, nhy]
    
    if type == 3    % Hyperbolic transform (velocity domain)
        nv = length(Param.v);           % Number of velocity parameters
        m = zeros(nt, nv);              % Initialize output model array
    else            % Linear or Parabolic transform (slowness domain)
        npx = length(Param.px);         % Number of x-slowness parameters
        npy = length(Param.py);         % Number of y-slowness parameters
        m = zeros(nt, npx, npy);        % Initialize output model array
    end
end

% =========================================================================
% RADON TRANSFORM IMPLEMENTATION
% =========================================================================
% The Radon transform is implemented using nested loops over all dimensions.
% For each model parameter (tau, px, py) or (tau, v), we compute the
% corresponding time shift and accumulate the contribution to the data.

switch type
    case 1  % LINEAR RADON TRANSFORM (Planar wavefronts)
        % Transform equation: t = tau + px*hx + py*hy
        % This models planar events in the data with constant slowness
        px = Param.px;  % x-slowness parameters
        py = Param.py;  % y-slowness parameters
        
        % Loop over all dimensions: intercept time, slowness parameters, and offsets
        for itau = 1:nt           % Loop over intercept time (model domain)
            for ipx = 1:length(px) % Loop over x-slowness parameters
                for ipy = 1:length(py) % Loop over y-slowness parameters
                    for ihx = 1:nhx    % Loop over x-offsets
                        for ihy = 1:nhy % Loop over y-offsets
                            % Compute travel time for this parameter combination
                            t = (itau-1)*dt + px(ipx)*hx(ihx) + py(ipy)*hy(ihy);
                            it = floor(t/dt) + 1;  % Convert to sample index
                            
                            % Check if computed time is within valid range
                            if it >= 1 && it <= nt
                                if operator == 1
                                    % Forward: Add model value to data at computed time
                                    d(it, ihx, ihy) = d(it, ihx, ihy) + m(itau, ipx, ipy);
                                else
                                    % Adjoint: Add data value to model at intercept time
                                    m(itau, ipx, ipy) = m(itau, ipx, ipy) + d(it, ihx, ihy);
                                end
                            end
                        end
                    end
                end
            end
        end

    case 2  % PARABOLIC RADON TRANSFORM (Parabolic wavefronts)
        % Transform equation: t = tau + px*hx^2 + py*hy^2
        % This models events with parabolic moveout (common in NMO-corrected data)
        px = Param.px;  % x-curvature parameters
        py = Param.py;  % y-curvature parameters
        
        for itau = 1:nt           % Loop over intercept time
            for ipx = 1:length(px) % Loop over x-curvature parameters
                for ipy = 1:length(py) % Loop over y-curvature parameters
                    for ihx = 1:nhx    % Loop over x-offsets
                        for ihy = 1:nhy % Loop over y-offsets
                            % Compute travel time with quadratic dependence on offsets
                            t = (itau-1)*dt + px(ipx)*hx(ihx)^2 + py(ipy)*hy(ihy)^2;
                            it = floor(t/dt) + 1;  % Convert to sample index
                            
                            if it >= 1 && it <= nt
                                if operator == 1
                                    % Forward: Accumulate model contribution to data
                                    d(it, ihx, ihy) = d(it, ihx, ihy) + m(itau, ipx, ipy);
                                else
                                    % Adjoint: Accumulate data contribution to model
                                    m(itau, ipx, ipy) = m(itau, ipx, ipy) + d(it, ihx, ihy);
                                end
                            end
                        end
                    end
                end
            end
        end

    case 3  % HYPERBOLIC RADON TRANSFORM (Hyperbolic wavefronts)
        % Transform equation: t = sqrt(tau^2 + (hx^2 + hy^2)/v^2)
        % This models events with hyperbolic moveout (common in pre-stack data)
        v = Param.v;  % Velocity parameters
        
        for itau = 1:nt           % Loop over zero-offset time
            for iv = 1:length(v)   % Loop over velocity parameters
                for ihx = 1:nhx    % Loop over x-offsets
                    for ihy = 1:nhy % Loop over y-offsets
                        % Compute distance from source to receiver
                        r2 = hx(ihx)^2 + hy(ihy)^2;  % Squared offset
                        
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
% Return the appropriate array based on operator direction
if operator == 1
    out = d;  % Forward: Return data space
else
    out = m;  % Adjoint: Return model space
end
end
