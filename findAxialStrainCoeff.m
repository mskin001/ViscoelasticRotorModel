function [E0,E1,C1,C2] = findAxialStrainCoeff(sigb)
global mat rim w

[Q, kappa, fi] = findMatPropConsts(1,1);

syms c1 c2 ez r
u = -mat.rho{1}*w^2*fi(1)*r^3 + c1*r^kappa + c2*r^-kappa + fi(2)*ez*r;
eT = u/r;
eR = diff(u,r);

SR = Q(1,3)*ez + Q(2,3)*eT + Q(3,3)*eR;
%% Find constants from boundaries
srInner = subs(SR,r,rim(1));
sInner = srInner - sigb(1);
C1 = solve(sInner,c1);

srOuter = subs(SR,{r,c1},[rim(end),C1]);
sOuter = srOuter - sigb(2);
C2 = solve(sOuter,c2);
C1 = subs(C1,c2,C2);

%% Find axial strain constants
sigZ = Q(1,1)*ez + Q(1,2)*eT + Q(1,3)*eR;
tempF = sigZ*r;
f = int(tempF,r,rim(1),rim(end));
f = subs(f,{c1,c2},[C1,C2]);
Ez = solve(f,ez);
Ez = double(Ez);

C1 = double(subs(C1,ez,Ez));
C2 = double(subs(C2,ez,Ez));

%% Check
InnerStress = subs(srInner,{c1,c2,ez},[C1,C2,Ez]) - sigb(1)
OuterStress = subs(srOuter,{c1,c2,ez},[C1,C2,Ez]) - sigb(2)


