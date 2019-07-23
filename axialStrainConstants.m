function [e0,e1] = axialStrainConstants(sigb, b)

global mat rim w;

Fw = 0;
F0 = Fw;
F1 = Fw;
Fc1 = Fw;
Fc2 = Fw;
Mw = Fw;
M1 = Fw;
M0 = Fw;
Mc1 = Fw;
Mc2 = Fw;
fw = zeros(2,1);
f1 = fw;
f0 = fw;
mw = fw;
m1 = fw;
m0 = fw;

for k = 1:length(rim)-1
  [~, kappa, fi] = findMatPropConsts(b,1);

  fw = - 0.25 * fi(9) * mat.rho{b} * w^2 * [-rim(k)^4; rim(k+1)^4]; 
  f1 = (1/3) * fi(12) * [-rim(k)^3; rim(k+1)^3];
  f0 = 0.5 * fi(13) * [-rim(k)^2; rim(k+1)^2];
  fc1 = fi(10)/(kappa+1) * [-rim(k)^(kappa+1); rim(k+1)^(kappa+1)];
  fc2 = fi(11)/(-kappa+1) * [-rim(k)^(-kappa+1); rim(k+1)^(-kappa+1)];
  
  Fw = Fw + sum(fw');
  F1 = F1 + sum(f1');
  F0 = F0 + sum(f0');
  Fc1 = Fc1 + sum(fc1');
  Fc2 = Fc2 + sum(fc2');
  
  mw = - (1/5) * fi(9) * mat.rho{b} * w^2 * [-rim(k)^5; rim(k+1)^5];
  m1 = 0.24 * fi(12) * [-rim(k)^4; rim(k+1)^4];
  m0 = (1/3) * fi(13) * [-rim(k)^3; rim(k+1)^3];
  mc1 = fi(10)/(kappa+2) * [-rim(k)^(kappa+2); rim(k+1)^(kappa+2)];
  mc2 = fi(11)/(-kappa+2) * [-rim(k)^(-kappa+2); rim(k+1)^(-kappa+2)];
  
  Mw = Mw + sum(mw');
  M0 = M0 + sum(m0');
  M1 = M1 + sum(m1');
  Mc1 = Mc1 + sum(mc1');
  Mc2 = Mc2 + sum(mc2');
  
end

syms si so ri ro c1 c2 ki ko f6i f6o f7i f7o sb1 sb2 e1 e0

inner = si*ri^2 + c1*ri^(ki-1) + c2*ri^(-ki-1) + f6i*e1*ri + f7i*e0 - sb1;
C1 = solve(inner,c1);

outer = so*ro^2 + C1*ro^(ko-1) + c2*ro^(-ko-1) + f6o*e1*ro + f7o*e0 - sb2;
C2 = solve(outer,c2);

C1 = subs(C1,c2,C2);

C2 = simplify(C2);
C1 = simplify(C1);

[~, ki, fi] = findMatPropConsts(b,1);
[~, ko, fo] = findMatPropConsts(b,(length(rim)-1));
ri = rim(1);
ro = rim(end);
sb1 = sigb(1);
sb2 = sigb(2);
si = mat.rho{1}*w^2*fi(6);
so = mat.rho{end}*w^2*fo(6);
f6i = fi(7);
f7i = fi(8);
f6o = fo(7);
f7o = fo(8);

C1 = subs(C1);
C2 = subs(C2);

Fz = Fw + F1 + F0 + Fc1*C1 + Fc2*C2;
Mz = Mw + M1 + M0 + Mc1*C1 + Mc2*C2;

e1 = solve(Fz,e1);
Mz = subs(Mz);
e0 = solve(Mz,e0);
e1 = subs(e1);

e0 = double(e0)
e1 = double(e1)





