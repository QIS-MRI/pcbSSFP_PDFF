% Function name: findLocalMinima
%
% Description: Find the local minima of the VARPRO residual at each
% voxel, and return them in a 3D array.
%
%Author: Diego Hernando
%Date created: March 18, 2008
%--------------------------------------------------------------------------

function [masksignal,resLocalMinima,numMinimaPerVoxel] = findLocalMinima( residual, threshold, masksignal )

L = size(residual,1);
sx = size(residual,2);
sy = size(residual,3);

dres = diff(residual,1,1);

if nargin < 3
  sumres = sqrt(squeeze(sum(residual,1)));
  sumres = sumres/max(max(sumres));
  masksignal = sumres>threshold;
end

resLocalMinima = zeros(1,sx,sy);
numMinimaPerVoxel = zeros(sx,sy,1);
for kx=1:sx
  for ky=1:sy
    if masksignal(kx,ky) > 0
      minres = min(residual(:,kx,ky));
      maxres = max(residual(:,kx,ky));

      temp = [0;squeeze(dres(:,kx,ky))];
      temp = temp<0 & circshift(temp,-1)>0 & residual(:,kx,ky)<minres+0.3*(maxres-minres);
      
      resLocalMinima(1:sum(temp),kx,ky) = find(temp);
      numMinimaPerVoxel(kx,ky) = sum(temp);
    end
  end
end
