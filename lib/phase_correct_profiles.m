function profiles = phase_correct_profiles(profiles, dataset, hWait)
%PHASE_CORRECT_PROFILES  Rotate to the real axis, then dictionary phase correction.
%
%   Repeats dataset.phasecorr.n_iters times (original scripts used 2).

if nargin < 3
    hWait = [];
end

tr = dataset.scan.tr;
te = dataset.scan.te;
fa = dataset.scan.fa;
pc_step = dataset.scan.pc_step;
rtrarr = dataset.phasecorr.rtrarr;
t2arr = dataset.phasecorr.t2arr;

freq_step_num = round(360 / pc_step);
[b0, rtr, t2] = ndgrid(linspace(0, 1e3/tr, freq_step_num), rtrarr, t2arr);
D1 = zeros(round(360/pc_step), numel(rtr));
for k = 1:numel(rtr)
    p = bSSFPAnalytic(rtr(k)*t2(k), t2(k), te, tr, fa, 0:pc_step:359, b0(k));
    D1(:, k) = conj(squeeze(p));
end

n_iters = dataset.phasecorr.n_iters;
for pass = 1:n_iters
    if ~isempty(hWait)
        sparcq_waitbar(hWait, (pass-1)/n_iters, ...
            sprintf('Phase correction  pass %d / %d', pass, n_iters));
    end
    profiles = profiles .* exp(-1i * angle(mean(profiles, 3)));
    [profiles, ~] = PhaseCorrection(profiles, conj(D1));
end
if ~isempty(hWait)
    sparcq_waitbar(hWait, 1, 'Phase correction  done');
end
end
