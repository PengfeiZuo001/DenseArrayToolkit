function config = loadConfig()
% loadConfig - Load all configuration parameters used by the program.
%
% Outputs:
%   config - A structure containing various configuration sub-structures.

    %% 1. Preprocessing parameters
    config.PreprocessingParam = struct();
    config.PreprocessingParam.tstart          = -60;   % start time (sec)
    config.PreprocessingParam.tend            = 600;   % end time (sec)
    config.PreprocessingParam.sig_leader      = 30;    % time before P
    config.PreprocessingParam.record_len      = 120;   % record length after P
    config.PreprocessingParam.lows            = 0.1;   % bandpass low corner (Hz)
    config.PreprocessingParam.highs           = 2.0;   % bandpass high corner (Hz)
    config.PreprocessingParam.resample_period = 0.1;   % resample period (sec)

    %% 2. Migration imaging parameters
    config.MigParam = struct();
    config.MigParam.is_ssa   = 0;     % SSA method or not
    config.MigParam.flow     = 0.1;   % min frequency (Hz)
    config.MigParam.fhigh    = 1.2;   % max frequency (Hz)
    config.MigParam.rank_p   = 8;     % rank
    config.MigParam.alpha    = 0.9;   % weight parameter
    config.MigParam.n_iter   = 20;    % number of iterations
    config.MigParam.gauss    = 2.5;   % Gaussian parameter
    config.MigParam.ph       = 5;     % phase parameter

    %% 3. Deconvolution parameters
    config.DeconvParam = struct();
    config.DeconvParam.gauss       = 5.0;    % Gaussian parameter
    config.DeconvParam.waterlevel  = 0.01;   % water-level parameter
    config.DeconvParam.itmax       = 100;    % max iteration
    config.DeconvParam.minderr     = 1e-5;   % minimum error
    config.DeconvParam.phaseshift  = 5;      % phase shift
    config.DeconvParam.verbose     = false;  % verbose output
    config.DeconvParam.radonfilter = false;  % use Radon filter or not

    %% 4. Radon Transform parameters
    config.RadonParam = struct();
    config.RadonParam.lows       = 0.1;      % bandpass low corner (Hz)
    config.RadonParam.highs      = 1.2;      % bandpass high corner (Hz)
    config.RadonParam.pmax       = 0.05;     % maximum slowness (s/km)
    config.RadonParam.pmin       = -0.05;    % minimum slowness (s/km)
    config.RadonParam.minTraces  = 60;       % minimum number of traces per event
    config.RadonParam.N1         = 30;       % CG iteration count
    config.RadonParam.N2         = 1;        % outer loop count
    config.RadonParam.plotRadon  = false;    % whether to plot Radon results

    %% 5. Array & event filtering parameters
    config.max_angle_diff  = 15; % max azimuth difference (deg)
    config.profile_length  = 4;  % profile length (km or other unit)

    %% 6. Global parameters
    config.dataFolder      = './data/event_waveforms_QBII'; % data folder path
    config.outputFolder    = './results';                 % output folder
    config.saveResults     = true;                        % whether to save results
    config.visualizeResults= true;                        % whether to visualize results
end