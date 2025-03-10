function idx = argmax(dat, dim)

if nargin == 1
[~,idx] = max(dat);
else
   [~,idx] = max(dat, [], dim); 
end