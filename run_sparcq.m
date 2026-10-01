%% SPARCQ  —  run on a prepared dataset
%
% Before this script:
%   1. Use prepare_data.m (this folder) to write a .mat with struct `dataset`
%      into example data/
%   2. Set prepared_name below to that filename
%
% Then run this file. A GUI asks for a background mask (Continue when done).
% Field-map residuals can take a while. Figures show ΔB0 and PDFF.
%
% Edit only the lines in USER SETTINGS.

setup_paths

%% USER SETTINGS
prepared_name = 'Liver_patients_ironOL_Pat2.mat';   % in example data/
save_results  = false;
show_figures  = true;

%% Load from example data/ (do not edit)
data_dir = example_data_dir();
prepared_file = fullfile(data_dir, prepared_name);

if ~isfile(prepared_file)
    error('File not found: %s\nRun prepare_data.m first (it writes into example data/).', prepared_file);
end

S = load(prepared_file);
if ~isfield(S, 'dataset')
    error('%s must contain a variable named dataset (from prepare_data.m).', prepared_file);
end
dataset = S.dataset;

fprintf('Loaded %s\n', prepared_file);
fprintf('  [%d x %d x %d]  TR=%.3f ms  FA=%.1f deg  type=%s\n', ...
    size(dataset.profiles,1), size(dataset.profiles,2), size(dataset.profiles,3), ...
    dataset.scan.tr, dataset.scan.fa, dataset.scan.type);

%% Mask  (GUI — adjust slider, then Continue)
mask = draw_intensity_mask(dataset.profiles);

%% Pipeline  (all other parameters come from dataset)
results = sparcq_pipeline(dataset, mask);

%% Figures
if show_figures
    mask_n = double(mask);
    mask_n(mask_n < 1) = nan;

    figure('Name', 'SPARCQ field map');
    imagesc(results.fm .* mask_n, [min(dataset.fieldmap.b0range), max(dataset.fieldmap.b0range)]);
    axis image off; colormap(gca, 'hot'); colorbar;
    title('\DeltaB_0 (Hz)');

    figure('Name', 'SPARCQ PDFF');
    imagesc(100 * results.FF .* mask_n, [0 100]);
    axis image off; colormap(gca, 'parula'); colorbar;
    title('PDFF (%)');
end

%% Optional: vial ROIs vs ground truth
% Left-click vials in the same order as gt_ff. Other mouse button to finish.
% ROIs are not saved.
plot_vials = true;
if plot_vials
    vial_radius = 5;   % pixels, passed to draw_vial_masks

    % Ground truth PDFF (%). Click vials in this order.
    % Paper names: commercial phantom = Calimetrix, custom phantom = Butters.
    commercial_phantom = [0, 2.5, 4.7, 7.3, 9.8, 14.0, 19.7, 29.4, 39.6, ...
                          49.9, 59.6, 68.9, 73.8, 79.1, 84.3, 89.0, 94.3, 100];
    custom_phantom     = [100, 100, 78, 49.9, 39.7, 30.8, 19.6, 14.8, 8.7, ...
                          6.4, 5.4, 3.2, 0, 0];

    gt_ff = commercial_phantom;   % or: gt_ff = custom_phantom;

    ff_pct = 100 * results.FF;
    ff_show = ff_pct;
    ff_show(double(mask) < 1) = nan;
    vial_masks = draw_vial_masks(ff_show, vial_radius);
    vial_stats = plot_vial_vs_gt(ff_pct, gt_ff, vial_masks);
    fprintf('Fit: measured = %.3f * GT + %.3f,  R^2 = %.3f\n', ...
        vial_stats.slope, vial_stats.intercept, vial_stats.R2);
end

%% Optional save
if save_results
    [p, n, ~] = fileparts(prepared_file);
    if isempty(p), p = pwd; end
    out_file = fullfile(p, [n '_results.mat']);
    save(out_file, 'results', 'dataset', 'mask', '-v7.3');
    fprintf('Saved %s\n', out_file);
end
