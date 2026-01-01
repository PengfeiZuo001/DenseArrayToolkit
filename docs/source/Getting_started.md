# Getting started

In the current version of DAT, Only SAC format files are surported and the program uses some information (e.g., KZTIME, DELTA, STLA, STLO, et al.) from the SAC headers. It is recommended to store all SAC files in one folder so they can be processed together. This section shows how to extract waveforms of earthquake from continuous recordings which are divided by days.

## Input files
We assume that the SAC files are all stored in the folder `DenseArrayToolkit-public/data/` and sac files are named as `yyyy.ddd.hh.mm.ss.xxx000.netcode.stacode...DPZ..SAC`, where `yyyy.ddd.hh.mm.ss.xxx000` is the origin time of the earthquake, `netcode` is the network code, `stacode` is the station code, and `DPZ` is vertical component.

Here is an example of the header of a SAC file (the file name is `2023.362.09.15.16.075.QB01.DPZ.SAC`):
```text
NPTS    = 30001
B       = 4.173980e+02             # necessary
E       = 1.017398e+03             # necessary
IFTYPE  = TIME SERIES FILE
LEVEN   = TRUE
DELTA   = 2.000000e-02             # necessary
DEPMIN  = -1.493024e-03
DEPMAX  = 1.382337e-03
DEPMEN  = 1.071187e-08
OMARKER = 0                        # necessary
T0MARKER= 477.4                    # necessary
KZDATE  = DEC 28 (362), 2023       # necessary
KZTIME  = 09:15:16.075             # necessary
KSTNM   = QB01                     # necessary
STLA    = 3.752070e+01             # necessary
STLO    = 9.138400e+01             # necessary
KEVNM   = -12345  -12345
EVLA    = 4.459070e+01             # necessary
EVLO    = 1.490094e+02             # necessary
EVDP    = 3.050100e+01             # necessary
DIST    = 4.805728e+03
AZ      = 2.814044e+02
BAZ     = 6.174802e+01
GCARC   = 4.323896e+01
LOVROK  = TRUE
USER0   = 4.632432e+02
NVHDR   = 6
LPSPOL  = FALSE
LCALDA  = TRUE
KCMPNM  = DPZ                       # necessary
KNETWK  = QB                        # necessary
MAG     = 6.500000e+00
```
We use ``read_SAC()`` to read the SAC file and get the header information. For detailed information, please refer to ``read_SAC()`` function.

## Parameters configuration

All default parameters are described in ``loadConfig.m`` file. Below is a brief description of the parameters. You can modify the parameters to suit your needs.

``` text
%%--------- Preprocessing parameters -------
config.PreprocessingParam = struct();
config.PreprocessingParam.tstart          = -60;   % start time (sec) relative to event
config.PreprocessingParam.tend            = 600;   % end time (sec) relative to event
config.PreprocessingParam.sig_leader      = 30;    % time before P-wave arrival (sec)
config.PreprocessingParam.record_len      = 120;   % record length after P-wave arrival (sec)
config.PreprocessingParam.lows            = 0.1;   % bandpass filter low corner frequency (Hz)
config.PreprocessingParam.highs           = 2.0;   % bandpass filter high corner frequency (Hz)
config.PreprocessingParam.resample_period = 0.1;   % resample period (sec)

%%--------- Deconvolution parameters -------
config.DeconvParam = struct();
config.DeconvParam.gauss       = 5.0;    % Gaussian parameter for deconvolution
config.DeconvParam.waterlevel  = 0.01;   % water-level parameter for deconvolution
config.DeconvParam.itmax       = 100;    % maximum number of iterations for deconvolution
config.DeconvParam.minderr     = 1e-5;   % minimum error for deconvolution
config.DeconvParam.phaseshift  = 5;      % phase shift for deconvolution
config.DeconvParam.verbose     = true;  % verbose output (true/false)
config.DeconvParam.radonfilter = false;  % use Radon filter (true/false)

%%--------- Radon Transform parameters -------
config.RadonParam = struct();
config.RadonParam.lows       = 0.1;      % bandpass filter low corner frequency (Hz)
config.RadonParam.highs      = 1.2;      % bandpass filter high corner frequency (Hz)
config.RadonParam.pmax       = 0.05;     % maximum slowness (s/km)
config.RadonParam.pmin       = -0.05;    % minimum slowness (s/km)
config.RadonParam.minTraces  = 60;       % minimum number of traces per event
config.RadonParam.N1         = 30;       % number of CG iterations
config.RadonParam.N2         = 1;        % number of outer loop iterations
config.RadonParam.plotRadon  = true;    % plot Radon results (true/false)

%%--------- Rank Reduction parameters (Off-the-grid reconstruction) -------
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

%%--------- Migration imaging parameters -------
config.MigParam = struct();
config.MigParam.is_ssa      = 0;     % use SSA method (1) or not (0)
config.MigParam.flow        = 0.1;   % minimum frequency for migration (Hz)
config.MigParam.fhigh       = 1.2;   % maximum frequency for migration (Hz)
config.MigParam.rank_p      = 8;     % rank parameter for SSA
config.MigParam.alpha       = 0.9;   % weight parameter for SSA
config.MigParam.n_iter      = 5;    % number of iterations for migration
config.MigParam.gauss       = 2.5;   % Gaussian parameter for migration
config.MigParam.phaseshift  = 5;     % phase parameter for migration
config.MigParam.plotMig     = false; % if true, plot migration results
%% CCP imaging parameters
config.CCPParam.imagingType = '3D';         % use 2D or 3D imaging, 2D imaging project the result onto the principal axis
config.CCPParam.plotCCP = false;            % if true, plot CCP results
config.CCPParam.smoothLength = 3;           % if greater than 0, apply smooth to CCP image

%%--------- Array & event filtering parameters -------
config.max_angle_diff  = 15; % max azimuth difference (deg)
config.profile_length  = 4;  % profile length (degree)

%%--------- Global parameters -------
config.dataFolder      = '../data/event_waveforms_BY'; % data folder path
config.outputFolder    = '../results';                 % output folder
config.saveResults     = true;                        % whether to save results
config.visualizeResults= true;                        % whether to visualize results
```