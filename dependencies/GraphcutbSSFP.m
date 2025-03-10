%% Function name: fw_i2cm1i_3pluspoint_hernando
%%
%% Description: Fat-water separation using regularized fieldmap formulation and graph cut solution. 
%%
%% Hernando D, Kellman P, Haldar JP, Liang ZP. Robust water/fat separation in the presence of large 
%% field inhomogeneities using a graph cut algorithm. Magn Reson Med. 2010 Jan;63(1):79-90.
%% 
%% Some properties:
%%   - Image-space
%%   - 2 species (water-fat)
%%   - Complex-fitting
%%   - Multi-peak fat (pre-calibrated)
%%   - Single-R2*
%%   - Independent water/fat phase
%%   - Requires 3+ echoes at arbitrary echo times (some choices are much better than others! see NSA...)
%%
%% Input: structures imDataParams and algoParams
%%   - imDataParams.images: acquired images, array of size[nx,ny,1,ncoils,nTE]
%%   - imDataParams.TE: echo times (in seconds)
%%   - imDataParams.FieldStrength: (in Tesla)
%%
%%   - algoParams.species(ii).name = name of species ii (string)
%%   - algoParams.species(ii).frequency = frequency shift in ppm of each peak within species ii
%%   - algoParams.species(ii).relAmps = relative amplitude (sum normalized to 1) of each peak within species ii
%%   Example
%%      - algoParams.species(1).name = 'water' % Water
%%      - algoParams.species(1).frequency = [0] 
%%      - algoParams.species(1).relAmps = [1]   
%%      - algoParams.species(2).name = 'fat' % Fat
%%      - algoParams.species(2).frequency = [3.80, 3.40, 2.60, 1.94, 0.39, -0.60]
%%      - algoParams.species(2).relAmps = [0.087 0.693 0.128 0.004 0.039 0.048]
%% 
%%   - algoParams.size_clique = 1; % Size of MRF neighborhood (1 uses an 8-neighborhood, common in 2D)
%%   - algoParams.range_r2star = [0 0]; % Range of R2* values
%%   - algoParams.NUM_R2STARS = 1; % Numbre of R2* values for quantization
%%   - algoParams.range_fm = [-400 400]; % Range of field map values
%%   - algoParams.NUM_FMS = 301; % Number of field map values to discretize
%%   - algoParams.NUM_ITERS = 40; % Number of graph cut iterations
%%   - algoParams.SUBSAMPLE = 2; % Spatial subsampling for field map estimation (for speed)
%%   - algoParams.DO_OT = 1; % 0,1 flag to enable optimization transfer descent (final stage of field map estimation)
%%   - algoParams.LMAP_POWER = 2; % Spatially-varying regularization (2 gives ~ uniformn resolution)
%%   - algoParams.lambda = 0.05; % Regularization parameter
%%   - algoParams.LMAP_EXTRA = 0.05; % More smoothing for low-signal regions
%%   - algoParams.TRY_PERIODIC_RESIDUAL = 0; % Take advantage of periodic residual if uniform TEs (will change range_fm)  
%%   - algoParams.residual: in case we pre-computed the fit residual (mostly for testing) 
%%
%% Output: structure outParams
%%   - outParams.species(ii).name: name of the species (taken from algoParams)
%%   - outParams.species(ii).amps: estimated water/fat images, size [nx,ny,ncoils] 
%%   - outParams.r2starmap: R2* map (in s^{-1}, size [nx,ny])
%%   - outParams.fieldmap: field map (in Hz, size [nx,ny])
%%
%%
%% Author: Diego Hernando
%% Date created: August 5, 2011
%% Date last modified: November 10, 2011
%% Modified for bSSFP by Berk Can ACikgoz
%% Date last modified: March 1, 2025

function fm = graphcut_bssfp( imDataParams, algoParams, fmguess )


DEBUG = 0;


% Check validity of params, and set default algorithm parameters if not provided
% [validParams,algoParams] = checkParamsAndSetDefaults_graphcut( imDataParams,algoParams );
% if validParams==0
%   disp(['Exiting -- data not processed']);
%   outParams = [];
%   return;
% end

% Get data dimensions
[sx,sy,N] = size(imDataParams.images);
% If more than one slice, pick central slice

% If precession is clockwise (positive fat frequency) simply conjugate data
if imDataParams.PrecessionIsClockwise <= 0 
  imDataParams.images = conj(imDataParams.images);
  imDataParams.PrecessionIsClockwise = 1;
end

% Check spatial subsampling option (speedup ~ quadratic SUBSAMPLE parameter)
% SUBSAMPLE = algoParams.SUBSAMPLE;
% if SUBSAMPLE > 1
%   images0 = imDataParams.images;
%   START = round(SUBSAMPLE/2);
%   [sx,sy] = size(images0(:,:,1,1,1));
%   allX = 1:sx;
%   allY = 1:sy;
%   subX = START:SUBSAMPLE:sx;
%   subY = START:SUBSAMPLE:sy;
%   imDataParams.images = images0(subX,subY,:,:,:);
% end

% Regularization parameter
lambda=algoParams.lambda;

% Spatially-varying regularization.  The LMAP_POWER applies to the
% sqrt of the curvature of the residual, and LMAP_POWER=2 yields
% approximately uniform resolution.
LMAP_POWER = algoParams.LMAP_POWER;
  
  
% LMAP_EXTRA: Extra flexibility for including prior knowledge into
% regularization. For instance, it can be used to add more smoothing
% to noise regions (by adding, eg a constant LMAP_EXTRA), or even to
% add spatially-varying smoothing as a function of distance to
% isocenter...
LMAP_EXTRA = algoParams.LMAP_EXTRA;

% Finish off with some optimization transfer -- to remove discretization
DO_OT = algoParams.DO_OT;
residual = algoParams.residual;

%save tempres.mat residual


% Setup the estimation, get the lambdamap,...
fms = linspace(algoParams.range_fm(1),algoParams.range_fm(2),algoParams.NUM_FMS);
dfm = fms(2)-fms(1);
lmap = getQuadraticApprox( residual, dfm );  
lmap = (sqrt(lmap)).^LMAP_POWER;
lmap = lmap + mean(lmap(:))*LMAP_EXTRA;

% Initialize the field map indices
cur_ind = ceil(length(fms)/2)*ones(size(imDataParams.images(:,:,1,1,1)));

% This is the core of the algorithm
fm = graphCutIterations(imDataParams,algoParams,residual,lmap,cur_ind );

% If we have subsampled (for speed), let's interpolate the field map
end





