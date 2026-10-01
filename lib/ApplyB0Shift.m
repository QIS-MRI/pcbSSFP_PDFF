function shifted_profile = ApplyB0Shift(profile, b0, tr)
%APPLYB0SHIFT  Circular shift of a phase-cycled profile by ΔB0 (Hz).

if size(profile, 1) < size(profile, 2)
    profile = profile.';
end
bw = 1e3 / tr;
b0 = b0 / bw;
Npc = length(profile);
pc_step = 360 / Npc;
pc = 0:pc_step:359;
[F, modenums] = FTMat(Npc-1, pc);
v = @(x) exp((1i) * 2 * pi * (2*modenums + 1) * x / 2);
shifted_profile = F' * diag(v(b0)) * F * profile;
end
