function [e0, e1] = axialStrainConstants(b,k)

global mat rim w;

syms e1 e2;

[Q, kappa, fi] = findMatPropConsts(b,k);

eta0 = -(Q(2,1) + 3*Q(2,3))*fi(1);
eta1 = (Q(2,1) + kappa*Q(2,3))*fi(2);
eta2 = (Q(2,1) - kappa*Q(2,3))*fi(3);
eta3 = (Q(2,1) + 2*Q(2,3))*fi(4) + Q(2,2);
eta4 = (Q(2,1) + Q(2,3))*fi(5) + Q(2,2);

for a = 1:length(rim) - 1
    fw = (1/4) * mat.rho * w^2 * eta0 * [-rim(a); rim(a+1)].^4;
    f1 = (1/3) * eta3 * e1 * [-rim(a); rim(a+1)].^3;
    f0 = (1/2) * eta4 * e0 * [-rim(a); rim(a+1)].^2;

    f = fw + f1 + f0 + (1/(kappa+1))*eta1*[-rim(a); rim(a+1)].^(kappa+1)...
        + (1/(kappa+1))*eta2*[-rim(a); rim(a+1)].^(-kappa+1);
    F(a:a+1) = F(a:a+1) + f;
    
    mw = (1/5) * mat.rho * w^2 * eta0 * [-rim(a); rim(a+1)].^5;
    m1 = (1/4) * eta3 * e1 * [-rim(a); rim(a+1)].^4;
    m0 = (1/3) * eta4 * e0 * [-rim(a); rim(a+1)].^3;

    m = mw + m1 + m0 + (1/(kappa+2))*eta1*[-rim(a); rim(a+1)].^(kappa+2)...
        - (1/(kappa+2))*eta2*[-rim(a); rim(a+1)].^(-kappa+2);
    M(a:a+1) = M(a:a+1)
end

