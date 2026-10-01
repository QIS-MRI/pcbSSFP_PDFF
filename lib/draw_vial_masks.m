function masks = draw_vial_masks(image, r)
%DRAW_VIAL_MASKS  Click vial centers; each click becomes a disk of radius r.
%
%   masks = draw_vial_masks(image, r)
%
%   Left-click each vial (same order as the ground-truth list).
%   Any other mouse button (or Enter) finishes.
%   masks is [Nx Ny nVials], 1 inside each circle.

if nargin < 2 || isempty(r)
    error('draw_vial_masks(image, r): give a circle radius in pixels.');
end

image = double(image);
limits = [0 100];
if all(~isfinite(image(:)))
    error('Image has no finite values to display.');
end
if ~any(image(:) >= 0 & image(:) <= 100)
    limits = [min(image(:), [], 'omitnan'), max(image(:), [], 'omitnan')];
end

hFig = figure('Name', 'Click vial centers', 'NumberTitle', 'off');
imagesc(image, limits);
axis image off
colormap(gca, 'parula')
colorbar
title(sprintf(['Left-click vial centers (radius %g px). ' ...
    'Other mouse button when done.'], r))
hold on

centers = zeros(0, 2);   % [row, col] for DrawCircle
while true
    [xi, yi, but] = ginput(1);
    if isempty(xi) || ~isequal(but, 1)
        break
    end
    row = round(yi);
    col = round(xi);
    centers(end+1, :) = [row, col]; %#ok<AGROW>
    plot(col, row, 'wo', 'MarkerSize', 8, 'LineWidth', 1.5)
    plot(col, row, 'kx', 'MarkerSize', 8, 'LineWidth', 1.2)
    text(col+3, row, sprintf('%d', size(centers,1)), ...
        'Color', 'w', 'FontWeight', 'bold', 'FontSize', 10)
end

if isempty(centers)
    close(hFig)
    error('No vials clicked.');
end

[Nx, Ny] = size(image);
nV = size(centers, 1);
masks = zeros(Nx, Ny, nV);
theta = linspace(0, 2*pi, 64);
for i = 1:nV
    masks(:, :, i) = DrawCircle([Nx, Ny], centers(i, :), r);
    plot(centers(i, 2) + r*cos(theta), centers(i, 1) + r*sin(theta), ...
        'w-', 'LineWidth', 1.2)
end
title(sprintf('%d vial(s), radius %g px', nV, r))
hold off
end
