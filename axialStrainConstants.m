function [e0, e1] = axialStrainConstants(sigb)

global mat rim w vari;

syms r eb ea cOne cTwo;
F = 0;
M = 0;

for b = 1:vari
  for k = 1:length(rim)-1
    [Q, kappa, fi] = findMatPropConsts(b,k);
    
    %Find C1 and C2 constants from radial boundary conditions
    sigR0 = -mat.rho{k}*w^2*fi(6)*rim(k)^2 + cOne*rim(k)^(kappa-1) + cTwo*rim(k)^(-kappa-1) + 2*fi(4)*eb*rim(k)...
      +fi(5)*ea - sigb(1);
    sigR1 = -mat.rho{k}*w^2*fi(6)*rim(k+1)^2 + cOne*rim(k+1)^(kappa-1) + cTwo*rim(k+1)^(-kappa-1) + 2*fi(4)*eb*rim(k+1)...
      +fi(5)*ea - sigb(2);
    C1 = solve(sigR0,cOne);
    sigR1 = subs(sigR1,cOne,C1);
    C2 = solve(sigR1,cTwo);
    C1 = subs(C1,cTwo,C2);
    
    ur = -mat.rho{k}*w^2*fi(1)*r^3 + C1*fi(2)*r^kappa + C2*fi(3)*r^-kappa + fi(4)*eb*r^2 + fi(5)*ea*r;

    eTheta = ur/r;
    er = diff(ur,r);

    sigZ = Q(2,1)*eTheta + Q(2,2)*eb*r + Q(2,2)*ea + Q(2,3)*er;
    zr = sigZ*r;
    m = sigZ*r^2;

    fz = int(zr,r,rim(k),rim(k+1));
    mz = int(m,r,rim(k),rim(k+1));

    F = F + fz;
    M = M + mz;
  end
end

e0 = solve(F,ea);
M = subs(M,ea,e0);
e1 = solve(M,eb);
e0 = subs(e0,eb,e1);
e1 = double(e1);
e0 = double(e0);
