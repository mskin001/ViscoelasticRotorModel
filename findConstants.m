function [C1,C2] = findConstants(sigb)

global mat rim w
syms c1 c2
b = 1; k = 1;
Q11 = mat.Q{b,k}(1,1);
Q13 = mat.Q{b,k}(1,3);
Q33 = mat.Q{b,k}(3,3);
kappa = sqrt(Q11/Q33);
fi3 = (3 * Q33 + Q13) / (Q33 * (9 - kappa^2));

sInner = -mat.rho{1}*w^2*fi3*rim(1)^2 + c1*rim(1)^(kappa-1) + c2*rim(1)^(-kappa-1) - sigb(1);
C1 = solve(sInner,c1);

sOuter = -mat.rho{1}*w^2*fi3*rim(end)^2 + C1*rim(end)^(kappa-1) + c2*rim(end)^(-kappa-1) - sigb(2);
C2 = solve(sOuter,c2);
C1 = subs(C1,c2,C2);

C1 = double(C1);
C2 = double(C2);