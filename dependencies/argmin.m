function idx = argmin(dat, dim)
if nargin==1
    [~,idx] = min(dat);
else
    [~, idx] = min(dat, [], dim);


end