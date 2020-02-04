t = linspace(1,3400,3400);
t0 = 0;

c1 = 33.3e-5;
c2 = 1.8391e-10;
c3 = 5.312e-7;
c4 = 1.6575e-6;
str = 5;

j = c1 + c2*(t-t0) + c3 * (1-exp((t0-t)*(c4/c3)));

elong = str * j;

h1 = 31.2;
matC1 = 4.11e-9;
matC2 = 1.82;
matD1 = -11.90;
matD2 = 1.6776;
matStr = [10.83 12.48 14.13 16.27 18.32 19.58];

for k = 1:length(matStr)
  h2 = log(1+h1) / (matD1 * matStr(k)^(matD2-1));
  eMat(k,1:length(t)) = (1/(h2*matStr(k))) * log(1 + (1+h1) * exp(h2*matC1*matStr(k)^(matC2+1)*t) - 1);
  S(k,1:length(t)) = eMat(k)/matStr(k);
  E(k,1:length(t)) = S(k,1:length(t)).^-1;
end