t = linspace(1,4000,4000);
t0 = zeros(1,length(t));

c1 = 3.3333e-5;
c2 = 1.8391e-10;
c3 = 5.3121e-7;
c4 = 1.6575e-6;
str = 79.25;

j = (c1 + c2*(t-t0)) + c3 * (1 - exp( (t0-t)*(c4/c3) ) );

elong = str * j;

h1 = 31.2;
matC1 = 4.11e-9;
matC2 = 1.82;
matD1 = -11.90;
matD2 = 1.6776;
matStr = 19.58;

for k = 1:length(matStr)
  h2 = (log(1+h1)) / (matD1 * matStr(k)^(matD2-1));
  eMat(k,1:length(t)) = (1/(h2*matStr(k))) * log(1 + (1+h1) * (exp(h2*matC1*t*matStr(k)^(matC2+1)) - 1));
    
  S(k,1:length(t)) = eMat(k)/matStr(k);
  E(k,1:length(t)) = S(k,1:length(t)).^-1;
  
  m = matC1*matStr(k)^matC2;
  q = matD1*matStr^matD2;
  e = m*t + q;
end


tao = [10 150 500];
psE1 = [113 0 1385 7784 9282];
psE2 = [401 0 3515 7072 10988];
psG = [58 315 5 258 636];

Elong = psE1(4) + psE1(1)*exp((-t/tao(1))) + psE1(2)*exp((-t/tao(2))) + psE1(3)*exp((-t/tao(3)));
Etrans = psE2(4) + psE2(1)*exp((-t/tao(1))) + psE2(2)*exp((-t/tao(2))) + psE2(3)*exp((-t/tao(3)));

transStrain = (ones(1,length(Etrans))*420)./Etrans;

figure(1), plot(t,Etrans)
figure(2), plot(t,Elong)
%-------------------------------------------------------------------------------
tao = [10 150 500];
ps = [113 0 1385; 401 0 3515];
pinf = [7784 7072];

E1 = pinf(1,1);
E2 = pinf(1,2);

for k = 1:length(tao)
  E1 = E1 + ps(1,k)*exp(-t/tao(k));
  E2 = E2 + ps(2,k)*exp(-t/tao(k));
end


E3 = E2;

S1 = 1/E1;
S2 = 1/E2;
S3 = 1/E3;

out = [S2,S3,S1];