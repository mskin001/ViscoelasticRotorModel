syms b a p r

A = (-p - 2*C)*a^2;
C = (a^2*p)/(2*(b^2 - a^2));

A = subs(A,C,C);
A = simplify(A);

sigRR = A/r^2 + 2*C;
sigRR = simplify(sigRR)
