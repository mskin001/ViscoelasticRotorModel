function [Q, kappa, fi] = findMatPropConsts(b,k)
% This function calculates the constants dependent on material properties.
% These are used for finding the stress and strain for each rim at at each
% velocity. The material property constants are dependent only on the
% stiffness of that rim. 
global mat
Q = mat.Q{b,k};

kappa = sqrt(Q(1,1)/Q(3,3));

% intermediate variables for calculating the stiffness matrix

fi(1) = 1/(Q(3,3) * (9 - kappa^2));
fi(2) = (Q(1,2)-Q(1,3))/(Q(3,3)-Q(2,2));
% fi(3) = 1/(Q(1,3) - (kappa * Q(3,3)));
% fi(4) = (Q(1,2) - 2*Q(2,3)) / (4*Q(3,3) - Q(1,1));
% fi(5) = (Q(1,2) - Q(2,3)) / (Q(3,3) - Q(1,1));
% fi(6) = fi(1) * (Q(1,3) + 3*Q(3,3));
% fi(7) = fi(4) * (Q(1,3) + 2*Q(3,3)) + Q(2,3);
% fi(8) = fi(5) * (Q(1,3) + Q(3,3)) + Q(2,3);
% fi(9) = fi(1) * (Q(1,2) + 3*Q(2,3));
% fi(10) = fi(2) * (Q(2,2) + kappa*Q(2,3));
% fi(11) = fi(3) * (Q(1,2) - kappa*Q(2,3));
% fi(12) = fi(4) * (Q(1,2) + 2*Q(2,3)) + Q(2,2);
% fi(13) = fi(5) * (Q(1,2) + Q(2,3)) + Q(2,2);

