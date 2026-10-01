function stats = plot_vial_vs_gt(ff_pct, gt_ff, masks)
%PLOT_VIAL_VS_GT  Vial-mean PDFF vs ground truth, with error bars and fit.
%
%   stats = plot_vial_vs_gt(ff_pct, gt_ff, masks)
%
%   ff_pct   PDFF map in percent
%   gt_ff    1 x nVials ground truth (same units, same click order)
%   masks    [Nx Ny nVials] from draw_vial_masks
%
%   Plot: GT on x, measured mean ± std on y, identity line, linear fit
%   on the vial means, and the equation y = ax + b.

gt_ff = gt_ff(:)';
nV = size(masks, 3);
if numel(gt_ff) ~= nV
    error('gt_ff has %d values but %d vials were clicked.', numel(gt_ff), nV);
end

ff_mean = zeros(1, nV);
ff_std  = zeros(1, nV);
for i = 1:nV
    vals = ff_pct(masks(:, :, i) > 0);
    vals = vals(isfinite(vals));
    if isempty(vals)
        error('Vial %d is empty (circle missed the map?).', i);
    end
    ff_mean(i) = mean(vals);
    ff_std(i)  = std(vals);
end

coeff = polyfit(gt_ff, ff_mean, 1);   % measured = a * GT + b
slope = coeff(1);
intercept = coeff(2);
yhat = polyval(coeff, gt_ff);
ss_res = sum((ff_mean - yhat).^2);
ss_tot = sum((ff_mean - mean(ff_mean)).^2);
if ss_tot == 0
    R2 = NaN;
else
    R2 = 1 - ss_res / ss_tot;
end

xline = linspace(min([0, gt_ff]), max([100, gt_ff]), 200);

figure('Name', 'Vial PDFF vs ground truth')
plot(xline, xline, 'k--', 'LineWidth', 1.5)
hold on
plot(xline, polyval(coeff, xline), 'r-', 'LineWidth', 1.5)
errorbar(gt_ff, ff_mean, ff_std, 'o', 'MarkerSize', 8, ...
    'LineWidth', 1.5, 'Color', [0 0.45 0.74], 'MarkerFaceColor', [0 0.45 0.74])
hold off
axis equal
grid on
xlim([-10 110])
ylim([-10 110])
xticks(0:20:100)
yticks(0:20:100)
xlabel('Ground truth PDFF (%)')
ylabel('Measured PDFF (%)')

if intercept >= 0
    eq = sprintf('y = %.3fx + %.3f', slope, intercept);
else
    eq = sprintf('y = %.3fx - %.3f', slope, abs(intercept));
end
if isfinite(R2)
    eq = sprintf('%s\nR^2 = %.3f', eq, R2);
end
text(5, 95, eq, 'FontSize', 12, 'FontName', 'Arial', ...
    'VerticalAlignment', 'top', 'BackgroundColor', 'w')
legend({'Identity', 'Fit', 'Mean \pm std'}, 'Location', 'southeast')
title('Vial PDFF vs ground truth')

stats.gt = gt_ff;
stats.mean = ff_mean;
stats.std = ff_std;
stats.slope = slope;
stats.intercept = intercept;
stats.R2 = R2;
end
