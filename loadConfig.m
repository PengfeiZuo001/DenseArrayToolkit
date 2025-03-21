function config = loadConfig()
% loadConfig - Load all configuration parameters used by the program.
%
% Outputs:
%   config - A structure containing various configuration sub-structures.

%% 1. Preprocessing parameters
config.PreprocessingParam = struct();
config.PreprocessingParam.tstart          = -60;   % start time (sec) relative to event
config.PreprocessingParam.tend            = 600;   % end time (sec) relative to event
config.PreprocessingParam.sig_leader      = 30;    % time before P-wave arrival (sec)
config.PreprocessingParam.record_len      = 120;   % record length after P-wave arrival (sec)
config.PreprocessingParam.lows            = 0.1;   % bandpass filter low corner frequency (Hz)
config.PreprocessingParam.highs           = 2.0;   % bandpass filter high corner frequency (Hz)
config.PreprocessingParam.resample_period = 0.1;   % resample period (sec)

%% 2. Deconvolution parameters
config.DeconvParam = struct();
config.DeconvParam.gauss       = 5.0;    % Gaussian parameter for deconvolution
config.DeconvParam.waterlevel  = 0.01;   % water-level parameter for deconvolution
config.DeconvParam.itmax       = 100;    % maximum number of iterations for deconvolution
config.DeconvParam.minderr     = 1e-5;   % minimum error for deconvolution
config.DeconvParam.phaseshift  = 5;      % phase shift for deconvolution
config.DeconvParam.verbose     = true;  % verbose output (true/false)
config.DeconvParam.radonfilter = false;  % use Radon filter (true/false)

%% 3. Radon Transform parameters
config.RadonParam = struct();
config.RadonParam.lows       = 0.1;      % bandpass filter low corner frequency (Hz)
config.RadonParam.highs      = 1.2;      % bandpass filter high corner frequency (Hz)
config.RadonParam.pmax       = 0.05;     % maximum slowness (s/km)
config.RadonParam.pmin       = -0.05;    % minimum slowness (s/km)
config.RadonParam.minTraces  = 60;       % minimum number of traces per event
config.RadonParam.N1         = 30;       % number of CG iterations
config.RadonParam.N2         = 1;        % number of outer loop iterations
config.RadonParam.plotRadon  = true;    % plot Radon results (true/false)

%% 4. Rank Reduction parameters (Off-the-grid reconstruction)
config.RankReductionParam.nx=11;   % number of grid points in x-direction
config.RankReductionParam.ny=14;   % number of grid points in y-direction
config.RankReductionParam.rank = 5;  % rank for rank reduction
config.RankReductionParam.niter = 5; % number of iterations
config.RankReductionParam.mode=1;    % mode for rank reduction
config.RankReductionParam.verb=true;    % verbosity flag
config.RankReductionParam.eps=0.00001; % epsilon for convergence
config.RankReductionParam.K=4;       % parameter K for rank reduction
config.RankReductionParam.flow=0.1;  % lower frequency bound
config.RankReductionParam.fhigh=1.2; % upper frequency bound
config.RankReductionParam.tmax=50;   % maximum time
config.RankReductionParam.plotRankReduction=false; % if true, plot rank reduction result
%% 5. Migration imaging parameters
config.MigParam = struct();
config.MigParam.is_ssa      = 0;     % use SSA method (1) or not (0)
config.MigParam.flow        = 0.1;   % minimum frequency for migration (Hz)
config.MigParam.fhigh       = 1.2;   % maximum frequency for migration (Hz)
config.MigParam.rank_p      = 8;     % rank parameter for SSA
config.MigParam.alpha       = 0.9;   % weight parameter for SSA
config.MigParam.n_iter      = 20;    % number of iterations for migration
config.MigParam.gauss       = 2.5;   % Gaussian parameter for migration
config.MigParam.phaseshift  = 5;     % phase parameter for migration
config.MigParam.plotMig     = false; % if true, plot migration results
%% CCP imaging parameters
config.CCPParam.imagingType = '3D';         % use 2D or 3D imaging, 2D imaging project the result onto the principal axis
config.CCPParam.plotCCP = false;            % if true, plot CCP results
config.CCPParam.smoothLength = 3;           % if greater than 0, apply smooth to CCP image

%% 5. Array & event filtering parameters
config.max_angle_diff  = 15; % max azimuth difference (deg)
config.profile_length  = 4;  % profile length (degree)

%% 6. Global parameters
config.dataFolder      = './data/event_waveforms_QBI'; % data folder path
config.outputFolder    = './results';                 % output folder
config.saveResults     = true;                        % whether to save results
config.visualizeResults= true;                        % whether to visualize results
end
