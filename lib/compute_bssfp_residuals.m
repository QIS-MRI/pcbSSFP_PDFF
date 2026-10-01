function res = compute_bssfp_residuals(profiles_lowres, mask, b0range, freqhz, amp, tr, fa, pc_step, hWait)
%COMPUTE_BSSFP_RESIDUALS  Low-resolution residual cube vs trial ΔB0.
%
%   res(x,y,ib0) = || measured - fitted 2-compartment bSSFP ellipse ||
%   Only voxels with mask>0 are fitted. Uses parfor over B0 (needs Parallel Toolbox).

if nargin < 9
    hWait = [];
end

profiles_lowres = profiles_lowres / abs(max(profiles_lowres(:)));
[Nx, Ny, ~] = size(profiles_lowres);
nB0 = numel(b0range);
res = zeros(Nx, Ny, nB0);

for x = 1:Nx
    for y = 1:Ny
        if mask(x, y) <= 0
            continue
        end
        p = squeeze(profiles_lowres(x, y, :));
        if angle(mean(p)) < 0
            p = -p;
        end
        parfor ib0 = 1:nB0
            b0 = b0range(ib0);
            [~, ~, ~, found_profile] = bSSFP_NonlinLS_Fit_WithPhase_Quant_MultComp( ...
                p, b0, b0 + freqhz, amp, pc_step, tr, fa);
            res(x, y, ib0) = norm(p - found_profile, 2);
        end
    end
    if ~isempty(hWait)
        sparcq_waitbar(hWait, x / Nx, sprintf('Field-map residuals  %.0f%%', 100 * x / Nx));
    end
end
end

function [Jm, F] = jf(x, J, f)
Jm = J(x);
F = f(x);
end

function e1 = getE1(a, b, alpha)
c = cosd(alpha);
e1 = (-b - c.*(b.*(a.^2)) + a + (a.*c)) ./ (a + (a.*c) - (b.*c) - (a.^2).*b);
end

function [sol, quant, res, found_profile] = bSSFP_NonlinLS_Fit_WithPhase_Quant_MultComp(p, phase, phase_fat, amps, pc_step, tr, fa)
% Local copy of the original residual fitter (algorithm unchanged).
if nargin < 5
    pc_step = 360 / length(p);
end
pc = (0:pc_step:359)';
pc = pc * pi / 180;
phase = phase * 2 * pi * (tr / 1e3);
pr = phase / 2;

phase_fat = phase_fat * 2 * pi * (tr / 1e3);
pr_fat = phase_fat * (1/2);

ellipse_canon = @(x) (x(3) * (((1 - x(1)*exp(1i*(-phase+pc)))) ./ (1 - x(2)*cos((-phase+pc))))) * exp(1i*pr);

u0 = [0.99; 0.7; 1];
u0_fat = [0.99; 0.7; 1];

phase_kern = [cos(pr)*eye(length(pc)), -sin(pr)*eye(length(pc)); ...
              sin(pr)*eye(length(pc)),  cos(pr)*eye(length(pc))];

j1c = @(x) x(3) * (exp(1i*(-phase+pc)) ./ (1 - x(2)*cos(-phase+pc)));
J1 = @(x) phase_kern * [real(j1c(x)); imag(j1c(x))];
j2c = @(x) -x(3) * ((1 - x(1)*exp(1i*(-phase+pc))) .* cos(-phase+pc) ./ ((1 - x(2)*cos(-phase+pc)).^2));
J2 = @(x) phase_kern * [real(j2c(x)); imag(j2c(x))];
j3c = @(x) -1 * (((1 - x(1)*exp(1i*(-phase+pc))) ./ (1 - x(2)*cos(-phase+pc))));
J3 = @(x) phase_kern * [real(j3c(x)); imag(j3c(x))];

jf_inline = cell(1, length(phase_fat));
for im = 1:length(phase_fat)
    phase_kern_fat = [cos(pr_fat(im))*eye(length(pc)), -sin(pr_fat(im))*eye(length(pc)); ...
                      sin(pr_fat(im))*eye(length(pc)),  cos(pr_fat(im))*eye(length(pc))];
    j1c_fat = @(x) x(6) * (exp(1i*(-phase_fat(im)+pc)) ./ (1 - x(5)*cos(-phase_fat(im)+pc)));
    j2c_fat = @(x) -x(6) * ((1 - x(4)*exp(1i*(-phase_fat(im)+pc))) .* cos(-phase_fat(im)+pc) ./ ...
        ((1 - x(5)*cos(-phase_fat(im)+pc)).^2));
    j3c_fat = @(x) -1 * (((1 - x(4)*exp(1i*(-phase_fat(im)+pc))) ./ (1 - x(5)*cos(-phase_fat(im)+pc))));
    J1_fat = @(x) phase_kern_fat * [real(j1c_fat(x)); imag(j1c_fat(x))];
    J2_fat = @(x) phase_kern_fat * [real(j2c_fat(x)); imag(j2c_fat(x))];
    J3_fat = @(x) phase_kern_fat * [real(j3c_fat(x)); imag(j3c_fat(x))];
    ellipse_canon_fat = @(x) (x(6) * (((1 - x(4)*exp(1i*(-phase_fat(im)+pc)))) ./ ...
        (1 - x(5)*cos((-phase_fat(im)+pc))))) * exp(1i*pr_fat(im));
    f = @(x) reimconcat(-ellipse_canon(x) - ellipse_canon_fat(x)) + reimconcat(p);
    J = @(x) [J1(x)'; J2(x)'; J3(x)'; J1_fat(x)'; J2_fat(x)'; J3_fat(x)']';
    jf_inline{im} = @(x) jf(x, J, f);
end

u0 = [u0; u0_fat];
u = u0;
options = optimoptions('lsqlin');
options.Display = 'off';
for it = 1:10
    Jm = zeros(length(reimconcat(p)), length(u));
    F = zeros(length(reimconcat(p)), 1);
    for im = 1:length(phase_fat)
        [Jmd, Fd] = jf_inline{im}(u);
        Jm = Jm + amps(im) * Jmd;
        F = F + amps(im) * Fd;
    end
    b = Jm * u - F;
    A = Jm;
    try
        u = lsqlin(A, b, ...
            [1 0 0 0 0 0; 0 1 0 0 0 0; 0 0 0 1 0 0; 0 0 0 0 1 0; ...
             -1 1 0 0 0 0; 0 0 0 -1 1 0; ...
             -1 0 0 0 0 0; 0 -1 0 0 0 0; 0 0 0 -1 0 0; 0 0 0 0 -1 0], ...
            [1 1 0.96 0.91 0 0 -0.82 -0.2 -0.82 -0.2], ...
            [], [], [], [], [], options);
    catch
        u = u0;
    end
end

sol = u;
found_profile = zeros(length(p), 1);
ellipse_canon = @(x) (x(3) * (((1 - x(1)*exp(1i*(-phase+pc)))) ./ (1 - x(2)*cos((-phase+pc))))) * exp(1i*pr);
for im = 1:length(amps)
    ellipse_canon_fat = @(x) (x(6) * (((1 - x(4)*exp(1i*(-phase_fat(im)+pc)))) ./ ...
        (1 - x(5)*cos((-phase_fat(im)+pc))))) * exp(1i*pr_fat(im));
    found_profile = found_profile + ((1/length(phase_fat)) * ellipse_canon(sol) + amps(im) * ellipse_canon_fat(sol));
end
f = @(x) reimconcat(-ellipse_canon(x) - ellipse_canon_fat(x)) + reimconcat(p);
res = norm(f(sol), 2);

a = abs(sol(1));
q = abs(sol(2));
M = abs(sol(3));
a_fat = abs(sol(4));
q_fat = abs(sol(5));
M_fat = abs(sol(6));
T2 = tr / log(1/a);
E1 = abs(getE1(a, q, fa));
T1 = abs(tr / log(1/E1));
pd = (M / ((1-E1)*sind(fa) / (1-E1*cosd(fa)-(sol(1)^2)*(E1-cosd(fa))))) / (exp(-(tr/2)/T2));
T2_fat = tr / log(1/a_fat);
E1_fat = abs(getE1(a_fat, q_fat, fa));
T1_fat = abs(tr / log(1/E1_fat));
pd_fat = (M_fat / ((1-E1_fat)*sind(fa) / (1-E1_fat*cosd(fa)-(sol(4)^2)*(E1_fat-cosd(fa))))) / (exp(-(tr/2)/T2_fat));
quant = [T1, T2, pd, T1_fat, T2_fat, pd_fat];
end
