syms b a p r x

C = (a^2*p)/(2*(b^2 - a^2));
A = (-p - 2*C)*a^2;

A = simplify(A);

sigRR = A/r^2 + 2*C;
sigRR = simplify(sigRR)

k = diag(5,5).*1
f = ones(5,1);

f = f + (x^2);

u = k\f;