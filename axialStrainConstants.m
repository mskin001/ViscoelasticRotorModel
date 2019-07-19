function [e0,e1] = axialStrainConstants(b)

global mat rim w;

B = zeros(3,3);
fw = zeros(2,1);
f1 = fw;
f0 = fw;
mw = fw;
m1 = fw;
m0 = fw;

[~, kone, fione] = findMatPropConsts(b,1);
[~, kend, fiend] = findMatPropConsts(b,length(rim)-1);

fw(1) = - 0.25 * fione(9) * mat.rho{1} * w^2 * -rim(1)^4;
fw(2) = - 0.25 * fiend(9) * mat.rho{end} * w^2 * rim(end)^4;

f1(1) = (1/3) * fione(12) * -rim(1)^3;
f1(2) = (1/3) * fiend(12) * rim(end)^3;

f0(1) = 0.5 * fione(13) * -rim(1)^2;
f0(2) = 0.5 * fiend(13) * rim(end)^2;

mw(1) = - (1/5) * fione(9) * mat.rho{1} * w^2 * -rim(1)^5;
mw(2) = - (1/5) * fiend(9) * mat.rho{1} * w^2 * rim(end)^5;

m1(1) = 0.24 * fione(12) * -rim(1)^4;
m1(2) = 0.24 * fiend(12) * rim(end)^4;

m0(1) = (1/3) * fione(13) * -rim(1)^3;
m0(2) = (1/3) * fiend(13) * rim(end)^3;

is = [-1 0; 0 1];
Gf = [rim(1)^(kone+1) rim(1)^(-kone+1); rim(end)^(kend+1) rim(end)^(-kend+1)];
Gm = [rim(1)^(kone+2) rim(1)^(-kone+2); rim(end)^(kend+2) rim(end)^(-kend+2)];
Kf = [1/(kone+1); 1/(-kend+1)];
Km = [1/(kone+2); 1/(-kend+2)];

A = is * Gf * (Kf * (Km.^-1)') * Gm^-1 * is^-1;

alpha = A*m1 - f1;
beta = A*m0 - f0;
omega = fw - A*mw;

e1 = (omega(2)*beta(1) - omega(1)*beta(2))/(beta(1)*alpha(2) - beta(2)*alpha(1));
e0 = omega(1)/beta(1) - (alpha(1)*e1)/beta(1);





