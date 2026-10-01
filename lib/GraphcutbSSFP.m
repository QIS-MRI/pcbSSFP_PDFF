%% Graph-cut field map from a precomputed bSSFP residual (Hernando 2010, adapted).

function fm = GraphcutbSSFP(imDataParams, algoParams, fmguess)

if imDataParams.PrecessionIsClockwise <= 0
    imDataParams.images = conj(imDataParams.images);
    imDataParams.PrecessionIsClockwise = 1;
end

LMAP_POWER = algoParams.LMAP_POWER;
LMAP_EXTRA = algoParams.LMAP_EXTRA;
residual = algoParams.residual;

fms = linspace(algoParams.range_fm(1), algoParams.range_fm(2), algoParams.NUM_FMS);
dfm = fms(2) - fms(1);
lmap = getQuadraticApprox(residual, dfm);
lmap = (sqrt(lmap)).^LMAP_POWER;
lmap = lmap + mean(lmap(:)) * LMAP_EXTRA;

cur_ind = ceil(length(fms)/2) * ones(size(imDataParams.images(:, :, 1, 1, 1)));
fm = graphCutIterations(imDataParams, algoParams, residual, lmap, cur_ind);
end
