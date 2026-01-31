function config = loadConfig()
% loadConfig - Load all configuration parameters used by the program.
%
% Outputs:
%   config - A structure containing various configuration sub-structures.

%% 1. Preprocessing parameters
config.PreprocessingParam                     = struct();
config.PreprocessingParam.tstart              = -60;        % start time (sec) relative to event
config.PreprocessingParam.tend                = 600;        % end time (sec) relative to event
config.PreprocessingParam.sig_leader          = 30;         % time before P-wave arrival (sec)
config.PreprocessingParam.record_len          = 120;        % record length after P-wave arrival (sec)
config.PreprocessingParam.lows                = 0.2;        % bandpass filter low corner frequency (Hz)
config.PreprocessingParam.highs               = 5.0;        % bandpass filter high corner frequency (Hz)
config.PreprocessingParam.resample_period     = 0.1;        % resample period (sec)

%% 2. Deconvolution parameters
config.DeconvParam                            = struct();
config.DeconvParam.gauss                      = 2.5;        % Gaussian parameter
config.DeconvParam.waterlevel                 = 0.01;       % water-level parameter
config.DeconvParam.itmax                      = 200;        % maximum number of iterations
config.DeconvParam.minderr                    = 1e-5;       % minimum error
config.DeconvParam.phaseshift                 = 5;          % phase shift
config.DeconvParam.verbose                    = true;       % verbose output
config.DeconvParam.radonfilter                = false;      % use Radon filter

%% 3. Radon Transform parameters
config.RadonParam                             = struct();
config.RadonParam.lows                        = 0.1;        % bandpass filter low corner frequency (Hz)
config.RadonParam.highs                       = 1.2;        % bandpass filter high corner frequency (Hz)
config.RadonParam.pmax                        = 0.05;       % maximum slowness (s/km)
config.RadonParam.pmin                        = -0.05;      % minimum slowness (s/km)
config.RadonParam.minTraces                   = 60;         % minimum number of traces
config.RadonParam.N1                          = 30;         % number of CG iterations
config.RadonParam.N2                          = 1;          % number of outer loop iterations
config.RadonParam.plotRadon                   = true;       % plot Radon results

%% 4. Rank Reduction parameters (Off-the-grid reconstruction)
config.RankReductionParam                     = struct();
config.RankReductionParam.nx                  = 11;         % number of grid points in x
config.RankReductionParam.ny                  = 14;         % number of grid points in y
config.RankReductionParam.rank                = 5;          % rank
config.RankReductionParam.niter               = 5;          % number of iterations
config.RankReductionParam.mode                = 1;          % mode
config.RankReductionParam.verb                = true;       % verbosity flag
config.RankReductionParam.eps                 = 1e-5;       % convergence threshold
config.RankReductionParam.K                   = 4;          % parameter K
config.RankReductionParam.flow                = 0.1;        % lower frequency bound
config.RankReductionParam.fhigh               = 1.2;        % upper frequency bound
config.RankReductionParam.tmax                = 50;         % maximum time
config.RankReductionParam.plotRankReduction   = false;      % plot rank reduction result

%% 5. Migration imaging parameters
config.MigParam                               = struct();
config.MigParam.bc                            = 1;          % boundary condition
config.MigParam.ssa                           = 0;          % use SSA method
config.MigParam.is_radon                      = 0;          % Radon on 3C waveforms

config.MigParam.flow                          = 0.1;        % minimum frequency (Hz)
config.MigParam.fhigh                         = 1.2;        % maximum frequency (Hz)
config.MigParam.itermax                       = 3;          % number of iterations 
config.MigParam.mu                            = 0.01;       % regularization parameter
config.MigParam.tol                           = 1e-6;       % tolerance
config.MigParam.gauss                         = config.DeconvParam.gauss;
config.MigParam.phaseshift                    = config.DeconvParam.phaseshift;
config.MigParam.src_type                      = 'p';        % source type
config.MigParam.fpeak                         = 1.2;        % peak frequency (unused)
config.MigParam.ispred                        = 0;          % predict waveforms
config.MigParam.t1                            = -5;         % start time of RF (sec)
config.MigParam.t2                            = 20;         % end time of RF (sec)
config.MigParam.minRatio                      = 0.6;        % minimum station ratio
config.MigParam.plotMig                       = false;      % plot migration results

%% 6. CCP imaging parameters
config.CCPParam                               = struct();
config.CCPParam.imagingType                   = '3D';       % 2D or 3D imaging
config.CCPParam.plotCCP                       = false;      % plot CCP results
config.CCPParam.smoothLength                  = 3;          % smoothing length

%% 7. Array & event filtering parameters
config.max_angle_diff                         = 15;         % max azimuth difference (deg)
config.profile_length                         = 4;          % profile length (degree)

%% 8. Global parameters
config.dataFolder                             = '../data/event_waveforms_QBI';  % data folder
config.outputFolder                           = './results';                     % output folder
config.saveResults                            = true;                           % save results
config.visualizeResults                       = true;                           % visualize results

end
