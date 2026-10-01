function [ff, fat_weights, water_weights] = Weights2FFFast(weights, b0, water_span)
%WEIGHTS2FFFAST  Fat fraction = weight mass outside a water-frequency window (Hz).

[Nx, Ny, ~] = size(weights);
b0arr = unique(b0);
ff = zeros(Nx, Ny);
fat_weights = zeros(Nx, Ny);
water_weights = zeros(Nx, Ny);
water_mask = (b0arr < (max(b0arr) - water_span)) + (b0arr > water_span) - 1;
water_comp = (b0arr > (max(b0arr) - water_span)) + (b0arr < water_span);

for x = 1:Nx
    for y = 1:Ny
        w = squeeze(weights(x, y, :));
        w = sum(reshape(abs(w), length(b0arr), length(b0)/length(b0arr)), 2);
        s = sum(abs(w));
        if s == 0
            continue
        end
        ff(x, y) = water_mask' * abs(w) / s;
        fat_weights(x, y) = water_mask' * abs(w);
        water_weights(x, y) = water_comp' * abs(w);
    end
end
end
