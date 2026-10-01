function fm = estimate_fieldmap(profiles_lowres, res, dataset)
%ESTIMATE_FIELDMAP  Hernando-style graph-cut on precomputed bSSFP residuals.

imDataParams = dataset.graphcut.imDataParams;
algoParams = dataset.graphcut.algoParams;

imDataParams.images = profiles_lowres;
algoParams.residual = permute(res, [3 1 2]);
algoParams.NUM_FMS = size(res, 3);
algoParams.range_fm = [min(dataset.fieldmap.b0range), max(dataset.fieldmap.b0range)];

[Nx, Ny, ~] = size(profiles_lowres);
fm = GraphcutbSSFP(imDataParams, algoParams, zeros(Nx, Ny));
end
