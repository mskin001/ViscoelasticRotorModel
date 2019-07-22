syms si so ri ro c1 c2 ki ko f6i f6o f7i f7o sb1 sb2 e1 e0

inner = si*ri^2 + c1*ri^(ki-1) + c2*ri^(-ki-1) + f6i*e1*ri + f7i*e0 - sb1;
C1 = solve(inner,c1);

outer = so*ro^2 + C1*ro^(ko-1) + c2*ro^(-ko-1) + f6o*e1*ro + f7o*e0 - sb2;
C2 = solve(outer,c2);

C1 = subs(C1,c2,C2);

C2 = simplify(C2);
C1 = simplify(C1);

sigr = si*ri^2 + C1*ri^(ki-1) + C2*ri^(-ki-1) + f6i*e1*ri + f7i*e0;
testInner = simplify(sigr)

sigr2 = so*ro^2 + C1*ro^(ko-1) + C2*ro^(-ko-1) + f6o*e1*ro + f7o*e0;
testOuter = simplify(sigr2)