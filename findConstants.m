function [C1, C2] = findConstants(sigb,b)

global mat w rim
% syms e1 e0 C1 C2
syms c1 c2 e1 e0

ri = rim(1);
ro = rim(end);

[~, ki, fi] = findMatPropConsts(b,1);
si = -mat.rho{1}*w^2*fi(6)*ri^2 + c1*ri^(ki-1) + c2*ri^(-ki-1) + fi(7)*e1*ri + fi(8)*e0 - sigb(1);
C1 = solve(si,c1);

[~, ko, fo] = findMatPropConsts(b,1);
so = -mat.rho{end}*w^2*fo(6)*ro^2 + C1*ro^(ko-1) + c2*ro^(-ko-1) + fo(7)*e1*ro + fo(8)*e0 - sigb(2);
C2 = solve(so,c2);

C1 = subs(C1,c2,C2);

e1 = 1;
e0 = 1;
sigr = subs(-mat.rho{1}*w^2*fi(6)*ri^2 + C1*ri^(ki-1) + C2*ri^(-ki-1) + fi(7)*e1*ri + fi(8)*e0);
testInner = double(subs(sigr))


sigr = subs(-mat.rho{end}*w^2*fo(6)*ro^2 + C1*ro^(ko-1) + C2*ro^(-ko-1) + fo(7)*e1*ro + fo(8)*e0);
testouter = double(subs(sigr))
