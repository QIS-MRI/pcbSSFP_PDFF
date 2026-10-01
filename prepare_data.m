%% SPARCQ data preparation
% Edit the USER SETTINGS block, then run this script.
% Raw profiles and the prepared dataset both live in example data/.
% Masking is not done here; use the GUI in run_sparcq.m.

clearvars -except profiles profiles_bias_corrected
close all

data_dir = example_data_dir();   % example data/

%% USER SETTINGS
raw_file  = 'PUT_YOUR_FILE.mat';                 % <-- input, in example data/
file_name = 'Commercial_phantom_example_dataset1.mat';   % <-- output, in example data/

% If the raw file uses another variable name (recon_cs, loaded_profiles, ...):
% profiles = squeeze(recon_cs(20,:,:,:));

%% Load from example data/
if ~exist('profiles','var')
    raw_path = fullfile(data_dir, raw_file);
    if ~isfile(raw_path)
        error('Raw data not found:\n  %s\nPut the .mat in example data/.', raw_path);
    end
    load(raw_path);
end

if ndims(profiles) ~= 3
    error('profiles must be [x y nPC]. Got size [%s].', num2str(size(profiles)));
end

%% Type defaults  (from the existing SPARCQ scripts)
% Pick one with `type` below. Override freq / amp / b0range after that if needed.

% Liver (in vivo)
typeDefaults.Liver.freq = [-3.800 -3.400 -2.600 -1.940 -0.390 0.600];
typeDefaults.Liver.amp  = [0.087 0.693 0.128 0.004 0.039 0.048];
typeDefaults.Liver.b0range = linspace(-400, 400, 101);   % Hz, ~8 Hz steps

% Knee (in vivo)
typeDefaults.Knee.freq = [-3.8 -3.4 -3.11 -2.7 -2.45 -1.93 -0.5 0.61];
typeDefaults.Knee.amp  = [0.087 0.568 0.058 0.092 0.058 0.027 0.038 0.073];
typeDefaults.Knee.b0range = linspace(-200, 200, 51);     % Hz, ~8 Hz steps

% Calimetrix phantom (commercial; 0.1096 ppm shift vs liver model)
typeDefaults.Calimetrix.freq = [-3.80, -3.40, -2.60, -1.94, -0.39, 0.60] - 0.1096;
typeDefaults.Calimetrix.amp  = [0.087 0.693 0.128 0.004 0.039 0.048];
typeDefaults.Calimetrix.b0range = linspace(-100, 100, 26);

% Butters phantom (custom-made)
typeDefaults.Butters.freq = [-3.950 -3.538 -3.266 -2.809 -2.595 -2.079 -0.751 -0.549 0.490];
typeDefaults.Butters.amp  = [0.09 0.607 0.047 0.058 0.036 0.025 0.031 0.038 0.069];
typeDefaults.Butters.b0range = linspace(-400, 400, 101);

%% USER SETTINGS  (change these per dataset)
% --- scan ---
tr = 3.4;                               % ms  (some datasets use 3.38)
te = tr/2;                              % ms
fa = 35;                                % degrees
field_strength = 2.89;                  % tesla
gyro = 42.6;                            % MHz/T  (scripts use 42.6, not 42.58)
precession_is_clockwise = 1;            % 1 = positive fat frequency

% --- fat spectrum / B0 search: 'Liver' | 'Knee' | 'Calimetrix' | 'Butters' ---
type = 'Calimetrix';
freq    = typeDefaults.(type).freq;     % ppm
amp     = typeDefaults.(type).amp;      % relative amplitudes
b0range = linspace(-100, 100, 26);  % Hz

% --- field-map downsample ---
scale_factor = 8;                       % 4 in recent scripts; older liver used 8

% --- phase correction ---
n_phase_corr_iters = 2;
phasecorr_rtrarr = 5;                   % T1/T2 for the 1-component dictionary
phasecorr_t2arr  = 80;                  % ms

% --- graph-cut (Hernando / GraphcutbSSFP) ---
% Defaults from main.m / mainEP_phantom.m / mainEP_knee.m
% (mainBCA.m used size_clique=2, lambda=1e-3, LMAP_EXTRA=1e-3, scale_factor=8)
gc_noise_bias_correction = 1;
gc_bssfp_flag = 1;
gc_size_clique = 3;                     % 1 = 8-neighborhood
gc_NUM_R2STARS = 10;
gc_NUM_ITERS = 100;
gc_SUBSAMPLE = 1;
gc_DO_OT = 0;
gc_LMAP_POWER = 2;
gc_lambda = 1e-2;
gc_LMAP_EXTRA = 1e-2;
gc_TRY_PERIODIC_RESIDUAL = 0;

% --- SPARCQ dictionary + NNLS ---
sparcq_rtrarr = 1:5:26;                 % T1/T2 samples
sparcq_t2arr = 80;                      % ms (single T2)
sparcq_freq_step_num = 36;              % off-resonance bins over 0..1/TR
sparcq_DF_res = 8;                      % Hz, kept for reference (unused in fit)
sparcq_reg_param = 0.15;                % Laplacian lambda (phantom and in vivo)

% --- SPARCQ fat-fraction from weights ---
water_span_wide = 120;                  % Hz  (Weights2FFFast)
water_span_narrow = 20;                 % Hz
ff_hi = 0.90;                           % overwrite SPARCQ with FFlookup above this
ff_lo = 0.02;                           % overwrite SPARCQ with FFlookup below this

% --- FFlookup relaxation (currently hard-coded inside FFlookup) ---
fflookup_t1w = 1e3;                     % ms
fflookup_t2w = 60;                      % ms
fflookup_t1f = 1e3;                     % ms
fflookup_t2f = 60;                      % ms
fflookup_nff = 101;                     % FF samples 0..1

%% Derived quantities  (do not edit unless you know why)
[Nx, Ny, nPC] = size(profiles);
pc_step = round(360/nPC);
pcnum = 360/pc_step;
pc_angles = 0:pc_step:359;

if ~isfield(typeDefaults, type)
    error('Unknown type ''%s''. Use Liver, Knee, Calimetrix, or Butters.', type);
end
if numel(freq) ~= numel(amp)
    error('freq and amp must have the same length.');
end

freqhz = freq * gyro * field_strength;  % ppm -> Hz
DF = linspace(0, 1000/tr, sparcq_freq_step_num);

%% Pack one structure
dataset.profiles = profiles;

dataset.scan.Nx = Nx;
dataset.scan.Ny = Ny;
dataset.scan.nPC = nPC;
dataset.scan.tr = tr;
dataset.scan.te = te;
dataset.scan.fa = fa;
dataset.scan.field_strength = field_strength;
dataset.scan.gyro = gyro;
dataset.scan.pc_step = pc_step;
dataset.scan.pcnum = pcnum;
dataset.scan.pc_angles = pc_angles;
dataset.scan.type = type;
dataset.scan.precession_is_clockwise = precession_is_clockwise;

dataset.fat.freq = freq;
dataset.fat.amp = amp;
dataset.fat.freqhz = freqhz;

dataset.phasecorr.n_iters = n_phase_corr_iters;
dataset.phasecorr.rtrarr = phasecorr_rtrarr;
dataset.phasecorr.t2arr = phasecorr_t2arr;

dataset.fieldmap.scale_factor = scale_factor;
dataset.fieldmap.b0range = b0range;

% Graph-cut: same fields GraphcutbSSFP / graphCutIterations read
dataset.graphcut.imDataParams.TR = tr * 1e-3;            % seconds
dataset.graphcut.imDataParams.FieldStrength = field_strength;
dataset.graphcut.imDataParams.PrecessionIsClockwise = precession_is_clockwise;

dataset.graphcut.algoParams.noise_bias_correction = gc_noise_bias_correction;
dataset.graphcut.algoParams.species(1).name = 'water';
dataset.graphcut.algoParams.species(1).frequency = 0;
dataset.graphcut.algoParams.species(1).relAmps = 1;
dataset.graphcut.algoParams.species(2).name = 'fat';
dataset.graphcut.algoParams.species(2).frequency = freq;
dataset.graphcut.algoParams.species(2).relAmps = amp;
dataset.graphcut.algoParams.bssfp_flag = gc_bssfp_flag;
dataset.graphcut.algoParams.size_clique = gc_size_clique;
dataset.graphcut.algoParams.NUM_R2STARS = gc_NUM_R2STARS;
dataset.graphcut.algoParams.range_fm = [min(b0range) max(b0range)];
dataset.graphcut.algoParams.NUM_FMS = numel(b0range);
dataset.graphcut.algoParams.NUM_ITERS = gc_NUM_ITERS;
dataset.graphcut.algoParams.SUBSAMPLE = gc_SUBSAMPLE;
dataset.graphcut.algoParams.DO_OT = gc_DO_OT;
dataset.graphcut.algoParams.LMAP_POWER = gc_LMAP_POWER;
dataset.graphcut.algoParams.lambda = gc_lambda;
dataset.graphcut.algoParams.LMAP_EXTRA = gc_LMAP_EXTRA;
dataset.graphcut.algoParams.TRY_PERIODIC_RESIDUAL = gc_TRY_PERIODIC_RESIDUAL;

dataset.sparcq.rtrarr = sparcq_rtrarr;
dataset.sparcq.t2arr = sparcq_t2arr;
dataset.sparcq.freq_step_num = sparcq_freq_step_num;
dataset.sparcq.DF_res = sparcq_DF_res;
dataset.sparcq.DF = DF;
dataset.sparcq.reg_param = sparcq_reg_param;
dataset.sparcq.water_span_wide = water_span_wide;
dataset.sparcq.water_span_narrow = water_span_narrow;
dataset.sparcq.ff_hi = ff_hi;
dataset.sparcq.ff_lo = ff_lo;

dataset.fflookup.t1w = fflookup_t1w;
dataset.fflookup.t2w = fflookup_t2w;
dataset.fflookup.t1f = fflookup_t1f;
dataset.fflookup.t2f = fflookup_t2f;
dataset.fflookup.nff = fflookup_nff;

dataset.typeDefaults = typeDefaults;

%% Save  (one variable in one .mat, into example data/)
save_path = fullfile(data_dir, file_name);
save(save_path, 'dataset', '-v7.3');

fprintf('Saved %s (variable: dataset)\n', save_path);
fprintf('  profiles [%d x %d x %d], TR=%.3f ms, FA=%.1f deg, B0=%.3f T, type=%s\n', ...
    Nx, Ny, nPC, tr, fa, field_strength, type);
fprintf('  b0range [%g, %g] Hz, %d steps; graph-cut lambda=%g, clique=%d\n', ...
    min(b0range), max(b0range), numel(b0range), ...
    dataset.graphcut.algoParams.lambda, dataset.graphcut.algoParams.size_clique);
fprintf('  SPARCQ lambda=%g, DF bins=%d, water spans %g / %g Hz\n', ...
    dataset.sparcq.reg_param, dataset.sparcq.freq_step_num, ...
    water_span_wide, water_span_narrow);
