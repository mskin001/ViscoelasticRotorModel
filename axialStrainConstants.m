function [e0] = axialStrainConstants(b,k)

global mat rim w;

syms r C1 C2 e0 e1;
F = 0;

[Q, kappa, fi] = findMatPropConsts(b,k);
  
ur = -mat.rho*w^2*fi(1)*r^3 + fi(2)*r^kappa + fi(3)*r^-kappa + fi(4)*eb*r^2 + fi(5)*ea*r;

eTheta = ur/r;
er = diff(ur,r);

sigZ = Q(2,1)*eTheta + Q(2,2)*ez + Q(2,3)*er;
zr = sigZ*r;
m = sigZ*r^2;

fz = int(zr,r,rim(k),rim(k+1));
mz = int(m,r,rim(k),rim(k+1));

F = F + fz;
M = M + mz;

e0 = solve(F,ea);
subs(M,ea,e0);
e1 = solve(M,eb);
subs(e0,eb,e1);
