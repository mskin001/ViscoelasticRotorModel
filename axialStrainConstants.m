function [e0] = axialStrainConstants(b,k)

global mat rim w;

syms r C1 C2 e0;

[Q, kappa, fi] = findMatPropConsts(b,k);

for a = 1:length(rim) - 1

  uz = -(mat.rho{a}*w^2*r^3)/((9-kappa^2)*Q(3,3)) + C1*fi(2)*r^kappa + C2*fi(3)*r^-kappa...
    + ((Q(2,2)-Q(3,1))/(Q(3,3)-Q(1,2)))*e0*r;
  
  eTheta = uz/r;
  er = diff(uz,r);
  
  sigZ = Q(2,1)*eTheta + Q(2,2)*e0 + Q(2,3)*er;
  sigZr = sigZ*r;
  
  fz = int(sigZr,r);
  
end
