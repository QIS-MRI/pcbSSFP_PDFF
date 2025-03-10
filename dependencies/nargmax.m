function [narg, idx] = nargmax(dat, n)

[sorted, sortidx] = sort(dat, "descend");
idx = sortidx(n);
narg = sorted(n);
end