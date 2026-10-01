function ff = FFlookup(cs, fm, dataset)
%FFLOOKUP  PDFF from chemical-shift vs field-map using a 2-compartment bSSFP lookup.

tr = dataset.scan.tr;
te = dataset.scan.te;
fa = dataset.scan.fa;
pc_step = dataset.scan.pc_step;
freqhz = dataset.fat.freqhz;
amp = dataset.fat.amp;
t1w = dataset.fflookup.t1w;
t2w = dataset.fflookup.t2w;
t1f = dataset.fflookup.t1f;
t2f = dataset.fflookup.t2f;
nff = dataset.fflookup.nff;

pcnum = 360 / pc_step;
f = zeros(pcnum, 1);
for k = 1:length(freqhz)
    f = f + amp(k) * bSSFPAnalytic(t1f, t2f, te, tr, fa, (0:pc_step:359)', freqhz(k));
end
w = bSSFPAnalytic(t1w, t2w, te, tr, fa, (0:pc_step:359)', 0);

ffarr = linspace(0, 1, nff);
offres_lookup = zeros(1, nff);
for k = 1:nff
    offres_lookup(k) = angle(mean((1-ffarr(k))*w + ffarr(k)*f)) * (1e3/tr) / pi;
end
ffarr = [ffarr, ffarr];
offres_lookup = [offres_lookup, 1e3/tr + offres_lookup];

offres_map = mod(cs - fm, 1e3/tr);
ff = zeros(size(cs));
for x = 1:size(cs, 1)
    for y = 1:size(cs, 2)
        ff(x, y) = ffarr(argmin(abs(offres_lookup - offres_map(x, y))));
    end
end
end
