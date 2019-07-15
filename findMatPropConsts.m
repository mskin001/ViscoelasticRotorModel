function [Q, kappa, fi] = findMatPropConsts(b,k)
% This function calculates the constants dependent on material properties.
% These are used for finding the stress and strain for each rim at at each
% velocity. The material property constants are dependent only on the
% stiffness of that rim. 
global mat

[r,c] = size(mat.Q{b,k});
Q = zeros(r,c);
for y = 1:r
    for x = 1:c
        Q(y,x) = mat.Q{b,k}(y,x);
    end
end

kappa = sqrt(Q(1,1)/Q(3,3));

% intermediate variables for calculating the stiffness matrix

fi(1) = 1/(Q(3,3) * (9 - kappa^2));
fi(2) = 1/(Q(1,3) + (kappa * Q(3,3)));
fi(3) = 1/(Q(1,3) - (kappa * Q(3,3)));
fi(4) = (Q(1,2) - 2*Q(2,3)) / (4*Q(3,3) - Q(1,1));
fi(5) = (Q(1,2) - Q(2,3)) / (Q(3,3) - Q(1,1));
fi(6) = (Q(1,3) + 3*Q(3,3))*fi(1);
fi(7) = fi(4) * (Q(1,3) + 2*Q(3,3)) + Q(2,3);
fi(8) = fi(5) * (Q(1,3) + Q(3,3)) + Q(2,3);

