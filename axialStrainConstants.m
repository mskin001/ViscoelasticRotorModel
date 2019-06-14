function [e0, e1] = axialStrainConstants(b,k, kappa)

global mat rim w;

Q11 = mat.Q{b,k}(1,1);
Q13 = mat.Q{b,k}(1,3);
Q12 = mat.Q{b,k}(1,2);
Q22 = mat.Q{b,k}(2,2);
Q23 = mat.Q{b,k}(2,3);
Q33 = mat.Q{b,k}(3,3);

fi0 = 1 / (Q33 * (9-kappa^2));
fi1 = 1 / (Q13 + kappa*Q33);
fi2 = 1 / (Q13 - kappa*Q33);
fi3 = (Q12 - 2*Q23) / (4*Q33 - Q11);
fi4 = (Q12 - Q23) / (Q33 - Q11);

fsig = (-3/3) * mat.rho{k} * w^2 * fi0 * (Q12 + Q23) * (-rim(k)^3 + rim(k+1)^3);
f1 = (1/2) * (fi3*(Q12 + 2*Q23) + Q22) * (-rim(k)^2 + rim(k+1)^2);
f0 = (1/1) * (fi4*(Q12 + Q23) + Q22) * (-rim(k)^1 + rim(k+1)^1);
fCT = (1/(kappa)) * (Q12 + kappa*Q23) * fi1 * (-rim(k)^(kappa) + rim(k+1)^(-kappa))...
      + (1/(-kappa)) * (Q12 + kappa*Q23) * fi2 * (-rim(k)^(kappa) + rim(k+1)^(-kappa));

msig = (-3/4) * mat.rho{k} * w^2 * fi0 * (Q12 + Q23)* (-rim(k)^4 + rim(k+1)^4);
m1 = (1/3) * (fi3*(Q12 + 2*Q23) + Q22) * (-rim(k)^3 + rim(k+1)^3);
m0 = (1/2) * (fi4*(Q12 + Q23) + Q22) * (-rim(k)^2 + rim(k+1)^2);
mCT = (1/(kappa+1)) * (Q12 + kappa*Q23) * fi1 * (-rim(k)^(kappa+1) + rim(k+1)^(-kappa+1))...
      + (1/(-kappa+1)) * (Q12 + kappa*Q23) * fi2 * (-rim(k)^(kappa+1) + rim(k+1)^(-kappa+1));

e0_inter = f1 / (m1*f0 - m0*f1);
e0 = e0_inter * (msig - (m1/f1) * fCT + mCT);

e1 = (f1^-1) * (-fsig - f0*e0 - fCT);
