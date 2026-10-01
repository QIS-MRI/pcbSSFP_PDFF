function circled = DrawCircle(sizeimg, c, r)
%DRAWCIRCLE  Disk mask of radius r centered at c = [row, col].

Nx = sizeimg(1);
Ny = sizeimg(2);
xc = c(1);
yc = c(2);
[jj, ii] = meshgrid(1:Ny, 1:Nx);
circled = double((ii - xc).^2 + (jj - yc).^2 < r^2);
end
