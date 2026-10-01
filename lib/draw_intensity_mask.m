function mask = draw_intensity_mask(profiles, fig_name)
%DRAW_INTENSITY_MASK  Interactive z-score mask. Click Continue when it looks right.
%
%   mask = draw_intensity_mask(profiles)
%
%   Slider is the z-score threshold used by IntensityMask (same as the original GUI).
%   The function blocks until Continue or the window is closed.

if nargin < 2 || isempty(fig_name)
    fig_name = 'Draw mask';
end

threshold = 0;
mask = IntensityMask(profiles, threshold);
mag = abs(sum(profiles, 3));

hFig = figure('Name', fig_name, 'NumberTitle', 'off', ...
    'MenuBar', 'none', 'ToolBar', 'none', ...
    'Position', [80 80 1000 560], 'Resize', 'off');
setappdata(hFig, 'mask', mask);
setappdata(hFig, 'threshold', threshold);

ax1 = axes('Parent', hFig, 'Position', [0.06 0.10 0.42 0.80]);
ax2 = axes('Parent', hFig, 'Position', [0.52 0.10 0.30 0.80]);
show_pair(ax1, ax2, mag, mask);

uicontrol(hFig, 'Style', 'text', 'String', 'Threshold (z-score)', ...
    'Position', [830 430 150 20], 'HorizontalAlignment', 'left');
sld = uicontrol(hFig, 'Style', 'slider', 'Min', -2, 'Max', 2, ...
    'Value', threshold, 'Position', [830 400 150 20]);
txt = uicontrol(hFig, 'Style', 'text', 'String', sprintf('%.3f', threshold), ...
    'Position', [830 370 150 20]);

addlistener(sld, 'Value', 'PostSet', @(~,~) on_slide());

uicontrol(hFig, 'Style', 'pushbutton', 'String', 'Continue', ...
    'Position', [830 320 150 36], 'FontWeight', 'bold', ...
    'Callback', @(~,~) uiresume(hFig));

uiwait(hFig);
if ~isvalid(hFig)
    error('Mask window was closed before Continue.');
end
mask = getappdata(hFig, 'mask');
close(hFig);

    function on_slide()
        threshold = get(sld, 'Value');
        set(txt, 'String', sprintf('%.3f', threshold));
        mask = IntensityMask(profiles, threshold);
        setappdata(hFig, 'mask', mask);
        setappdata(hFig, 'threshold', threshold);
        show_pair(ax1, ax2, mag, mask);
    end
end

function show_pair(ax1, ax2, mag, mask)
imagesc(ax1, mag .* mask); axis(ax1, 'image', 'off'); colormap(ax1, 'gray');
title(ax1, 'Masked magnitude');
imagesc(ax2, mask); axis(ax2, 'image', 'off'); colormap(ax2, 'gray');
title(ax2, 'Mask');
end
