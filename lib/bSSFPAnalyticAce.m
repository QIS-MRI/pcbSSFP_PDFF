function [signal, a, b, M, theta] = bSSFPAnalyticAce(t1, t2, te, tr, alpha, pc, off_res)
%BSSFPANALYTICACE  Same analytic bSSFP as bSSFPAnalytic, extra ellipse parameters.

e1 = exp(-tr/t1);
e2 = exp(-tr/t2);
alpha = alpha * pi / 180;
d = (1 - e1*cos(alpha) - (e2^2)*(e1 - cos(alpha)));
b = e2 * (1-e1) * (1+cos(alpha)) / d;
a = e2;
M = (1-e1) * sin(alpha) / d;
delt = pc * pi / 180;
theta0 = 2 * pi * off_res * (tr * 1e-3);
theta = theta0 - delt;
phi = theta0 * te / tr;
signal = (M * (1 - a*exp(1i*theta)) ./ (1 - b*cos(theta))) * exp(-1i*phi);
signal = conj(signal) * exp(-te/t2);
M = M * exp(-te/t2);
end
