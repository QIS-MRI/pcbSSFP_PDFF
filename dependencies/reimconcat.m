function concat_res = reimconcat(signal, dim)

    if nargin == 1
        signal_re = real(signal);
        signal_im = imag(signal);
        concat_res = [signal_re;signal_im];
    else
        signal_re = real(signal);
        signal_im = imag(signal);
        concat_res = cat(dim, signal_re, signal_im);
    end

end