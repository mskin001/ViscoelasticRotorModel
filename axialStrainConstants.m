function [e0, e1] = axialStrainConstants(sigR, b,k, kappa)

global mat rim w;

fw = (1/4) * mat.rho * w^2 * eta0 * r^4;
f1 = (1/3) * eta3 * e1 * r^3;
f0 = (1/2) * eta4 * e0 * r^2;

f = fw + f1 + f0 + (1/(k+1))*eta1*r^(k+1) + (1/(k+1))*eta2*r^(-k+1);

mw = (1/5) * mat.rho * w^2 * eta0 * r^5;
m1 = (1/4) * eta3 * e1 * r^4;
m0 = (1/3) * eta4 * e0 * r^3;

m = mw + m1 + m0 + (1/(k+2))*eta1*r^(k+2) - (1/(k+2))*eta2*r^(-k+2);




