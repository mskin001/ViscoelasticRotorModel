function [out] = Glass_Epoxy_Marty2019()

global t

tao = [10 150 500];
ps = [113 0 1385; 401 0 3515; 58 315 5];
pinf = [7784 7072 258];

E1 = pinf(1,1);
E2 = pinf(1,2);
G = pinf(1,3);

for k = 1:length(tao)
  E1 = E1 + ps(1,k)*exp(-t/tao(k));
  E2 = E2 + ps(2,k)*exp(-t/tao(k));
  G = G + ps(3,k)*exp(-t/tao(k));
end

E1 = E1*10^6;
E2 = E2*10^6;
G = G*10^6;

E3 = E2;

S1 = 1/E1;
S2 = 1/E2;
S3 = 1/E3;

out = [S1, S2, G];