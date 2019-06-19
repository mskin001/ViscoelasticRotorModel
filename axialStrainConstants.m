function [e0, e1] = axialStrainConstants(sigR, b,k, kappa)

global mat rim w;

Q11 = mat.Q{b,k}(1,1);
Q13 = mat.Q{b,k}(1,3);
Q12 = mat.Q{b,k}(1,2);
Q22 = mat.Q{b,k}(2,2);
Q23 = mat.Q{b,k}(2,3);
Q33 = mat.Q{b,k}(3,3);

fi0 = 1 / (Q33 * (9-kappa^2));
fi3 = (Q12 - 2*Q23) / (4*Q33 - Q11);
fi4 = (Q12 - Q23) / (Q33 - Q11);
fi5 = (Q13 + 3*Q33)*fi0;
fi6 = fi3 * (Q13 + 2*Q33) + Q23;
fi7 = fi4 * (Q13 + Q33) + Q23;

% e0 indipendent intermediate terms
a1 = -(1/2) * (Q22 - (fi7*(Q12+Q23)*kappa^-1)) * [-rim(k)^2; rim(k+1)^2];
a1 = a1(2) + a1(1); 
a2 = (1/2) * (Q12+Q23) * kappa^-1 * sigR * [-rim(k)^2; rim(k+1)^2];
a2 = a2(2) + a2(1);
a3 = (1/4) * mat.rho{b} * w^2 * (fi5*(Q12+3*Q23)*kappa^-1 + fi6*(Q12+3*Q23))*[-rim(k)^4; rim(k+1)^4];
a3 = a3(2) + a3(1);
a4 = (1/3) * (Q22 - fi6*(Q12-2*Q23)*kappa^-1) * [-rim(k)^3; rim(k+1)^3]; % * e1
a4 = a4(2) + a4(1);
a5 = (1/2) * fi4 * (Q12 + Q23) * [-rim(k)^2; rim(k+1)^2];
a5 = a5(2) + a5(1);
a6 = (1/3) * fi3 * (Q12 + 2*Q23) * [-rim(k)^3; rim(k+1)^3];
a5 = a5(2) + a5(1);

% e1 intermediate terms
b1 = (1/4) * (Q22 - fi6*(Q12 - 2*Q23)*kappa^-1) * [-rim(k)^4; rim(k+1)^4];
b1 = b1(2) + b1(1);
b2 = (1/3) * (Q12+Q23) * kappa^-1 * sigR *[-rim(k)^3; rim(k+1)^3];
b2 = b2(2) + b2(1);
b3 = (1/5) * mat.rho{b} * w^2 * (fi5*(Q12+3*Q23)*kappa^-1 + fi6*(Q12+3*Q23)) * [-rim(k)^5; rim(k+1)^5];
b3 = b3(2) + b3(1);
b4 = (1/3) * (Q22 - fi7*(Q12+Q23)*kappa^-1)*[-rim(k)^3; rim(k+1)^3]; % * e0
b4 = b4(2) + b4(1);
b5 = (1/3) * fi4 * (Q12+Q23)*[-rim(k)^3; rim(k+1)^3];
b5 = b5(2) + b5(1);
b6 = (1/4) * (Q12 + 2*Q23)*[-rim(k)^4; rim(k+1)^4];
b6 = b6(2) + b6(1);

e1 = (1/(1-(b4*a4)/(b1*a1))) * ((1/b1) * ((b2 + b3 + b5 + b6) + (b4/a1)*(a2 + a3 + a5 + a6)));
e0 = (1/a1) * (a2 + a3 + a4*e1 + a5 + a6);
