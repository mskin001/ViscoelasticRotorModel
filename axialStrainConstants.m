function [e0,e1] = axialStrainConstants(b,k,C1,C2)

global mat rim w;

Fw = zeros(2,1);
F0 = Fw;
F1 = Fw;
Mw = Fw;
M1 = Fw;
M0 = Fw;
fw = Fw;
f1 = Fw;
f0 = Fw;
mw = Fw;
m1 = Fw;
m0 = Fw;
for k = 1:length(rim)-1
  [~, kappa, fi] = findMatPropConsts(b,1);

  fw = - 0.25 * fi(9) * mat.rho{b} * w^2 * [-rim(k)^4; rim(k+1)^4]; 
  
  f1 = (1/3) * fi(12) * [-rim(k)^3; rim(k+1)^3];
  f0 = 0.5 * fi(13) * [-rim(k)^2; rim(k+1)^2];
  Fw = Fw + fw;
  F1 = F1 + f1;
  F0 = F0 + f0;
  
  mw = - (1/5) * fi(9) * mat.rho{b} * w^2 * [-rim(k)^5; rim(k+1)^5];
  m1 = 0.24 * fi(12) * [-rim(k)^4; rim(k+1)^4];
  m0 = (1/3) * fi(13) * [-rim(k)^3; rim(k+1)^3];
  Mw = Mw + mw;
  M0 = M0 + m0;
  M1 = M1 + m1;
end

is = [-1 0; 0 1];
Gf = [rim(1)^(kappa+1) rim(1)^(-kappa+1); rim(end)^(kappa+1) rim(end)^(-kappa+1)];
Gm = [rim(1)^(kappa+2) rim(1)^(-kappa+2); rim(end)^(kappa+2) rim(end)^(-kappa+2)];
Kf = [1/(kappa+1); 1/(-kappa+1)];
Km = [1/(kappa+2); 1/(-kappa+2)];

A = is * Gf * (Kf * (Km.^-1)') * Gm^-1 * is^-1;

alpha = A*m1 - f1;
beta = A*m0 - f0;
omega = fw - A*mw;

e1 = (omega(2)*beta(1) - omega(1)*beta(2))/(beta(1)*alpha(2) - beta(2)*alpha(1));
e0 = omega(1)/beta(1) - (alpha(1)*e1)/beta(1);





