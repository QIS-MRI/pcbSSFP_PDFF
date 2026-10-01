function [weights, FF] = sparcq_fit(profiles, fm, mask, dataset, hWait)
%SPARCQ_FIT  Dictionary NNLS per voxel, then spectral fat fraction.

if nargin < 5
    hWait = [];
end

tr = dataset.scan.tr;
te = dataset.scan.te;
fa = dataset.scan.fa;
pc_step = dataset.scan.pc_step;
rtrarr = dataset.sparcq.rtrarr;
t2arr = dataset.sparcq.t2arr;
DF = dataset.sparcq.DF;
reg_param = dataset.sparcq.reg_param;

[b0, rtr, t2] = ndgrid(DF, rtrarr, t2arr, 0);
parameters1 = [b0(:), rtr(:)];
D1 = zeros(round(360/pc_step), numel(rtr));
for k = 1:numel(rtr)
    p = bSSFPAnalyticAce(rtr(k)*t2(k), t2(k), te, tr, fa, (0:pc_step:359)', b0(k));
    D1(:, k) = squeeze(p);
end

profiles = profiles / max(abs(profiles(:)));
[Nx, Ny, ~] = size(profiles);
weights = zeros(Nx, Ny, size(D1, 2));

for x = 1:Nx
    for y = 1:Ny
        if mask(x, y) <= 0
            continue
        end
        p1 = squeeze(profiles(x, y, :));
        p1 = ApplyB0Shift(p1, -fm(x, y), tr);
        p1 = normalize(p1, 'norm');
        if angle(mean(p1)) * (1e3/tr) / pi < 0
            w = NNLS_Laplace(reimconcat(-p1), reg_param, reimconcat(D1), parameters1);
        else
            w = NNLS_Laplace(reimconcat(p1), reg_param, reimconcat(D1), parameters1);
        end
        weights(x, y, :) = abs(w);
    end
    if ~isempty(hWait)
        sparcq_waitbar(hWait, x / Nx, sprintf('SPARCQ matching  %.0f%%', 100 * x / Nx));
    end
end

FF1 = Weights2FFFast(abs(weights), squeeze(parameters1(:, 1)), dataset.sparcq.water_span_wide);
FF2 = Weights2FFFast(abs(weights), squeeze(parameters1(:, 1)), dataset.sparcq.water_span_narrow);
FF = 0.5 * (FF1 + FF2);
FF(~isfinite(FF)) = 0;
end
