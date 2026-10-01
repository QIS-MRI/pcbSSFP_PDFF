% Function name: createExpansionGraphVARPRO_fast
% Description: adjacency graph for the fieldmap update (min-cut == best 2^Q neighbor)
% Author: Diego Hernando, March 18, 2008

function [A,Aunsc] = createExpansionGraphVARPRO_fast( residual, dfm, lambda, size_clique, cur_ind, step )

sx = size(residual,2);
sy = size(residual,3);
L = size(residual,1);
s = sx*sy;
num_nodes = s + 2;
sA = [num_nodes, num_nodes];

offset = [0:(s-1)]'*L;
step_ind = cur_ind + step;

if strcmp(computer,'GLNX86')
  maxA = 1e5;
else
  maxA = 1e6;
end

valsh = zeros(1,s);
valsv = zeros(s,1);

factor = lambda*dfm^2;
x = 1:sx;
y = 1:sy;
[Y,X] = meshgrid(y,x);

allIndCross = [];
allValsCross = [];
for dx = -size_clique:size_clique
  for dy = -size_clique:size_clique
    dist = sqrt(dx^2+dy^2);
    if dist>0
      validmapi = X+dx>=1 & X+dx<=sx & Y+dy>=1 & Y+dy<=sy;
      validmapj = X-dx>=1 & X-dx<=sx & Y-dy>=1 & Y-dy<=sy;

      curfactor = min(factor(validmapi),factor(validmapj));

      a = curfactor.*(1./dist*(cur_ind(validmapi)-cur_ind(validmapj)).^2);
      b = curfactor.*(1./dist*(cur_ind(validmapi)-step_ind(validmapj)).^2);
      c = curfactor.*(1./dist*(step_ind(validmapi)-cur_ind(validmapj)).^2);
      d = curfactor.*(1./dist*(step_ind(validmapi)-step_ind(validmapj)).^2);

      temp = zeros(1,s);
      temp(validmapi) = max(0,c-a);
      valsh = valsh + temp;

      temp = zeros(s,1);
      temp(validmapi) = max(0,a-c);
      valsv = valsv + temp;

      temp = zeros(1,s);
      temp(validmapj) = max(0,d-c);
      valsh = valsh + temp;

      temp = zeros(s,1);
      temp(validmapj) = max(0,c-d);
      valsv = valsv + temp;

      S = 1:s;
      Sh = 1 + (S(validmapi(:)));
      Sv = 1 + (S(validmapj(:)));
      indcross = sub2ind(sA, Sh,Sv );
      temp = b+c-a-d;

      allIndCross = [allIndCross,indcross];
      allValsCross = [allValsCross,temp(:)'];
    end
  end
end

ind1 = sub2ind(sA, ones(1,num_nodes-2),2:(num_nodes-1) );
ind2 = sub2ind(sA,2:(num_nodes-1) , num_nodes + zeros(1,num_nodes-2));

temp0 = residual(cur_ind(:)+offset(:));
valid_ind = (step_ind(:)>=1 & step_ind(:)<=L);
temp1 = zeros(s,1);
temp1(valid_ind~=0) = residual(step_ind(valid_ind~=0)+offset(valid_ind~=0));
curmaxA = max(max(temp0),max([valsh.';valsv;allValsCross.']));
infty = curmaxA;
temp1(valid_ind==0) = infty;

indAll = [ind1 ind2 allIndCross];
valuesAll = [valsh + reshape(max(temp1-temp0,0),1,s), valsv.' + reshape(max(0,temp0-temp1),1,s), allValsCross];

[indAllSort,sortIndex] = sort(indAll);
valsort = valuesAll(sortIndex);
[xind,yind] = ind2sub([num_nodes num_nodes], indAllSort);
A = sparse(xind,yind,valsort,num_nodes,num_nodes,length(valsort));

if nargout>1
  Aunsc = A;
end

A = round(A*maxA/curmaxA);
A(A<0) = 0;
