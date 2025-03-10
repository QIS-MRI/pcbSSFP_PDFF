function [Fm, modenums] = FTMat(num_of_modes, cycles)
cycles = cycles*pi/180;

if mod(num_of_modes,2)
    modes = (1:1:num_of_modes) - ((num_of_modes+1)/2);
else
    modes = (1:1:num_of_modes) - (num_of_modes+1)/2;
end
modenums = modes;
% modes = -5:1:5;

Fm = zeros(length(modes), length(cycles));

for i = 1:length(modes)
    for j = 1:length(cycles)
        Fm(i,j) = exp(1i*modes(i)*cycles(j));
    end
end

Fm = Fm/sqrt(length(cycles));

end

