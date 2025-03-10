function [peak, loc] = getNthPeak(p, n, infsup)

if infsup == "inf"
    p = -1*p;
end

    [pks,locs] = findpeaks(p);
    try
    peak = pks(n);
    loc = locs(n);
    catch
        peak = [];
        loc = [];
    end
end
