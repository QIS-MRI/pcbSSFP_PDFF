function concat_res = reimconcat(signal, dim)
%REIMCONCAT  Stack real and imaginary parts (rows, or along dim).

if nargin == 1
    concat_res = [real(signal); imag(signal)];
else
    concat_res = cat(dim, real(signal), imag(signal));
end
end
