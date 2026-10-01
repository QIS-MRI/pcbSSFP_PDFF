function [profiles, phase_offset] = PhaseCorrection(profiles, dictionary)
%PHASECORRECTION  Voxel-wise phase offset from a real-valued dictionary match.

D = dictionary;
M = pinv(real(D' * D));
[Nx, Ny, ~] = size(profiles);
phase_offset = zeros(Nx, Ny);
corrected = zeros(size(profiles));

for x = 1:Nx
    for y = 1:Ny
        p = squeeze(profiles(x, y, :));
        h = D' * p;
        c = h.' * M * h;
        phase_offset(x, y) = angle(c);
    end
end

phase_offset(isnan(phase_offset)) = 0;

for x = 1:Nx
    for y = 1:Ny
        a = phase_offset(x, y) / 2;
        p = squeeze(profiles(x, y, :));
        corrected(x, y, :) = p * exp(-1i * a);
        phase_offset(x, y) = -a;
    end
end
profiles = corrected;
end
