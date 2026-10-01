function [weights, bestFitSignal, regMatrix] = NNLS_Laplace(profile, lambda, lib, parameters)
%NNLS_LAPLACE  Non-negative least squares with 2-D Laplacian on (off-resonance, T1/T2).

arrayX = unique(parameters(:, 1));
arrayY = unique(parameters(:, 2));

acquiredSignalWithReg = [profile; zeros(size(lib, 2), 1)];

sz = length(arrayX) * length(arrayY);
lgX = length(arrayX);

dfVector = ones(sz-1, 1);
dfVector(lgX:lgX:end) = 0;
dfVectorWrap = zeros(lgX * length(arrayY) - lgX + 1, 1);
dfVectorWrap(1:lgX:end) = 1;
rtrVector = ones(sz - lgX, 1);
diagVector = -4 * ones(sz, 1);
diagVector(1:lgX) = -3;
diagVector(end-lgX:end) = -3;

dfLaplacian = diag(diagVector, 0) + diag(dfVector, 1) + diag(dfVector, -1) ...
    + diag(dfVectorWrap, lgX-1) + diag(dfVectorWrap, -(lgX-1));
rtrLaplacian = diag(rtrVector, lgX) + diag(rtrVector, -lgX);

regMatrix = lambda * (dfLaplacian + rtrLaplacian);
dictionaryWithReg = [lib; regMatrix];

weights = lsqnonneg(dictionaryWithReg, acquiredSignalWithReg);
bestFitSignal = lib * weights;
end
