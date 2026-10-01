function [low, mask_low] = downsample_profiles(profiles, mask, scale_factor)
%DOWNSAMPLE_PROFILES  Spatial downsample of profiles and mask for field mapping.

nPC = size(profiles, 3);
sz = imresize(abs(profiles(:, :, 1)), 1/scale_factor);
low = zeros([size(sz), nPC], 'like', profiles);
for ipc = 1:nPC
    low(:, :, ipc) = imresize(squeeze(profiles(:, :, ipc)), 1/scale_factor);
end
mask_low = imresize(double(mask), 1/scale_factor, 'nearest') > 0;
end
