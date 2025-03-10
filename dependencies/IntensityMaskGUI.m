function intense_mask = IntensityMaskGUI(profiles)
    % Create a figure for the GUI
    hFig = figure('Position', [100, 100, 1000, 600], 'MenuBar', 'none', ...
                   'Name', 'Intensity Masking', 'NumberTitle', 'off', ...
                   'Resize', 'off');

    % Create axes for the masked image and the mask image
    ax1 = axes('Position', [0.1, 0.1, 0.4, 0.8]);  % Masked image
    ax2 = axes('Position', [0.55, 0.1, 0.4, 0.8]);  % Mask image
    threshold = 0;  % Initial threshold value

    % Initial mask and masked image
    intense_mask = IntensityMask(profiles, threshold);
    maskedImage = abs(sum(profiles, 3)) .* intense_mask;

    % Display the initial images
    imshow(maskedImage, [], 'Parent', ax1);
    title(ax1, 'Masked Image');
    imshow(intense_mask, [], 'Parent', ax2);
    title(ax2, 'Mask Image');

%     % Create a panel for the threshold slider and controls
%     uipanel('Position', [0.75, 0.1, 0.2, 0.8], 'Title', 'Threshold Adjustment');

    % Create a slider for threshold adjustment
    uicontrol('Style', 'text', 'Position', [600, 520, 150, 20], ...
              'String', 'Adjust Threshold:');
    thresholdSlider = uicontrol('Style', 'slider', 'Min', -2, 'Max', 2, ...
                                 'Value', threshold, 'Position', [600, 490, 150, 20]);
    thresholdValueText = uicontrol('Style', 'text', 'Position', [675, 460, 50, 20], ...
                                    'String', num2str(get(thresholdSlider, 'Value')));

    % Create a button to save the current mask
    uicontrol('Style', 'pushbutton', 'String', 'Save Mask', ...
              'Position', [600, 430, 150, 30], 'Callback', @saveMask);

    % Slider callback to update threshold value and apply mask automatically
    addlistener(thresholdSlider, 'Value', 'PostSet', @(src, event) updateImages(thresholdSlider, profiles, ax1, ax2, thresholdValueText));

    function updateImages(sliderHandle, profiles, ax1, ax2, textHandle)
        threshold = get(sliderHandle, 'Value');
        set(textHandle, 'String', num2str(threshold));  % Update the threshold value display
        intense_mask = IntensityMask(profiles, threshold);  % Update the mask based on the new threshold
        maskedImage = abs(sum(profiles, 3)) .* intense_mask;  % Update the masked image

        % Update the masked image and mask image displays
        axes(ax1);
        imshow(maskedImage, [], 'Parent', ax1);
        title(ax1, 'Masked Image');

        axes(ax2);
        imshow(intense_mask, [], 'Parent', ax2);
        title(ax2, 'Mask Image');
    end

    function saveMask(~, ~)
        % Save the current mask as a variable in the workspace
        assignin('base', 'intense_mask', intense_mask);
        assignin('base', 'threshold', threshold);
        msgbox('Mask saved as variable ''intense_mask'' in the workspace.','Success');
    end
end
