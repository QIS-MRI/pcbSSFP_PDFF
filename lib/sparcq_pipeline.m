function results = sparcq_pipeline(dataset, mask)
%SPARCQ_PIPELINE  Full PDFF pipeline using a prepare_data `dataset` struct.
%
%   results = sparcq_pipeline(dataset, mask)
%
%   dataset  struct from prepare_data.m (profiles + scan/fat/phasecorr/...)
%   mask     [x y] logical/double foreground (from draw_intensity_mask)
%
%   results.FF      PDFF map in [0, 1]
%   results.FF_GC   lookup PDFF (used at very low / very high FF)
%   results.fm      ΔB0 map (Hz), full resolution
%   results.weights SPARCQ NNLS weights [x y nDict]
%   results.cs      chemical-shift map (Hz)
%   results.mask    mask used

profiles = dataset.profiles;
tr = dataset.scan.tr;
te = dataset.scan.te;
fa = dataset.scan.fa;
pc_step = dataset.scan.pc_step;
scale = dataset.fieldmap.scale_factor;
b0range = dataset.fieldmap.b0range;
mask = double(mask > 0);

hWait = sparcq_waitbar([], 0, 'Phase correction...');
cleanup = onCleanup(@() sparcq_waitbar(hWait, 'close'));

%% 1. Phase correction
profiles = phase_correct_profiles(profiles, dataset, hWait);

%% 2. Downsample for field mapping
hWait = sparcq_waitbar(hWait, 0, 'Field-map residuals...');
profiles(isnan(profiles)) = 0;
[profiles_low, mask_low] = downsample_profiles(profiles, mask, scale);

%% 3. Voxel-wise bSSFP residuals vs ΔB0
res = compute_bssfp_residuals(profiles_low, mask_low, b0range, ...
    dataset.fat.freqhz, dataset.fat.amp, tr, fa, pc_step, hWait);
res(~isfinite(res)) = 0;
res(res < 0) = 0;

%% 4. Graph-cut field map, then upsample
hWait = sparcq_waitbar(hWait, 0.92, 'Graph-cut field map...');
fm_low = estimate_fieldmap(profiles_low, res, dataset);
fm = imresize(fm_low, scale);

%% 5. Chemical shift + lookup PDFF
hWait = sparcq_waitbar(hWait, 0.96, 'Lookup PDFF...');
cs = mod(angle(mean(profiles, 3)) * (1e3/tr) / pi, 1e3/tr);
FF_GC = FFlookup(cs, fm, dataset);

%% 6. SPARCQ dictionary matching
hWait = sparcq_waitbar(hWait, 0, 'SPARCQ matching...');
[weights, FF_sparcq] = sparcq_fit(profiles, fm, mask, dataset, hWait);

FF = FF_sparcq;
FF(FF > 1) = 1;
FF(FF_GC > dataset.sparcq.ff_hi) = FF_GC(FF_GC > dataset.sparcq.ff_hi);
FF(FF_GC < dataset.sparcq.ff_lo) = FF_GC(FF_GC < dataset.sparcq.ff_lo);
FF(mask < 1) = 0;

results.FF = FF;
results.FF_GC = FF_GC;
results.fm = fm;
results.cs = cs;
results.weights = weights;
results.mask = mask;
sparcq_waitbar(hWait, 1, 'Done.');
end
