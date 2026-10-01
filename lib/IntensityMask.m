function [intense_mask, intense_mask_nand] = IntensityMask(profiles, threshold, op)
%INTENSITYMASK  Foreground from z-scored |sum of phase cycles|.

if nargin < 3
    op = 10;
end

complex_sum = abs(sum(profiles, 3));
feature = standardize_image(complex_sum);
intense_mask = double(bwareaopen(feature > threshold, op));
intense_mask_nand = intense_mask;
intense_mask_nand(intense_mask == 0) = nan;
end

function data = standardize_image(data)
sz = size(data);
data = data(:);
data = data - mean(data);
data = data / std(data);
data = reshape(data, sz);
end
