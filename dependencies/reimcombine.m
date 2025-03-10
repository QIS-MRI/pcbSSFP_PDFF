function complex = reimcombine(x)
Nx = size(x,1); 
realx = x(1:Nx/2, :); 
imagx = x((Nx/2)+1:end, :);
complex = realx + 1i*imagx;
end