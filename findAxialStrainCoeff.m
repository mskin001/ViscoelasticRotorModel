function [E,C] = findAxialStrainCoeff(sigb)
global mat rim w vari






%% Find axial strain constants
for b = 1:vari
  syms c1 c2 ez r fi1 fi2 k Q13 Q23 Q33 
  u = -mat.rho{1}*w^2*fi1*r^3 + c1*r^k + c2*r^-k + fi2*ez*r;
  eT = u/r;
  eR = diff(u,r);

  SR = Q13*ez + Q23*eT + Q33*eR;
  %% Find constants from boundaries
  [Q, kappa, fi] = findMatPropConsts(b,1);
  r = rim(1);
  fi1 = fi(1);
  fi2 = fi(2);
  k = kappa;
  Q13 = Q(1,3);
  Q23 = Q(2,3);
  Q33 = Q(3,3);
  srInner = subs(SR) - sigb(1);
  C1 = solve(srInner,c1);

  [~, kappa, fi] = findMatPropConsts(b,length(rim)-1);
  r = rim(end);
  fi1 = fi(1);
  fi2 = fi(2);
  k = kappa;
  Q13 = Q(1,3);
  Q23 = Q(2,3);
  Q33 = Q(3,3);
  c1 = C1;
  srOuter = subs(SR) - sigb(2);
  C2 = solve(srOuter,c2);
  C1 = subs(C1,c2,C2);
  
  syms r
  [Q, kappa, fi] = findMatPropConsts(b,1);
  u = -mat.rho{1}*w^2*fi(1)*r^3 + C1*r^kappa + C2*r^-kappa + fi(2)*ez*r;
  eT = u/r;
  eR = diff(u,r);

  sigZ = Q(2,1)*ez + Q(2,2)*eT + Q(2,3)*eR;
  tempF = sigZ*r;
  f = int(tempF,r,rim(1),rim(2));
  Ez = solve(f,ez);
  E(b) = double(Ez);

  C1 = double(subs(C1,ez,E(b)));
  C2 = double(subs(C2,ez,E(b)));
  C(b,:) = [C1;C2];
  
  [Q, kappa, fi] = findMatPropConsts(b,1);
  r = rim(1);
  fi1 = fi(1);
  fi2 = fi(2);
  k = kappa;
  Q13 = Q(1,3);
  Q23 = Q(2,3);
  Q33 = Q(3,3);
  c1 = C(b,1);
  c2 = C(b,2);
  ez = E(b);
  srInner = double(subs(SR) - sigb(1))
end

%% Check
% [Q, kappa, fi] = findMatPropConsts(1,1);
% r = rim(1);
% fi1 = fi(1);
% fi2 = fi(2);
% k = kappa;
% Q13 = Q(1,3);
% Q23 = Q(2,3);
% Q33 = Q(3,3);
% c1 = C(1);
% c2 = C(2);
% ez = Ez;
% srInner = subs(SR) - sigb(1)
% 
% [Q, kappa, fi] = findMatPropConsts(1,1);
% r = rim(end);
% fi1 = fi(1);
% fi2 = fi(2);
% k = kappa;
% Q13 = Q(1,3);
% Q23 = Q(2,3);
% Q33 = Q(3,3);
% c1 = C(1);
% c2 = C(2);
% ez = Ez;
% srInner = subs(SR) - sigb(end)
