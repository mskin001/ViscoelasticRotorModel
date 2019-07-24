function [E0,E1,C1,C2] = findAxialStrainCoeff(sigb)
global mat rim w

[Q, kappa, fi] = findMatPropConsts(1,1);

syms c1 c2 e0 e1 r
%% Find constants from boundaries
sInner = -mat.rho{1}*w^2*fi(5)*rim(1)^2 + c1*rim(1)^(kappa-1) + c2*rim(1)^(-kappa-1)...
          + fi(7)*e1*rim(1) + fi(8)*e0 - sigb(1);
C1 = solve(sInner,c1);

sOuter = -mat.rho{1}*w^2*fi(5)*rim(end)^2 + C1*rim(end)^(kappa-1) + c2*rim(end)^(-kappa-1)...
          + fi(7)*e1*rim(end) + fi(8)*e0 - sigb(2);
C2 = solve(sOuter,c2);
C1 = subs(C1,c2,C2);

%% Find axial strain constants
u = -mat.rho{1}*w^2*fi(1)*r^3 + C1*fi(2)*r^kappa + C2*fi(3)*r^-kappa + fi(4)*e1*r^2 ...
      + fi(5)*e0;
eT = u/r;
eR = diff(u,r);
eZ = e0 + e1*r;

sigZ = Q(2,1)*eT + Q(2,2)*eZ + Q(2,3)*eR;
tempF = sigZ*r;
f = int(tempF,r,rim(1),rim(end));
E1 = solve(f,e1);
sigZ = subs(sigZ,e1,E1);
tempM = sigZ * r^2;
m = int(tempM,r,rim(1),rim(end));
E0= solve(m,e0);
E0 = double(E0)

E1 = subs(E1,e0,E0);
E1 = double(E1)

C1 = double(subs(C1,{e0,e1},[E0,E1]))
C2 = double(subs(C2,{e0,e1},[E0,E1]))
