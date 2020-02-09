t = linspace(1,1800,1800);
A = 5.43e-2;
n = 0.0117;
e0 = 0;
fm = e0 + A*t.^n;

Em = 1673.7e6;
num = 5.47e12;
Ek = 3719.6e6;
nuk = 4.36e11;
sig0 = 2e6;

% sig0 = 1;
% Em = 0.273;
% num = 0.219;
% Ek = 0.001;
% nuk = -0.276;
% t = [0,100,200,400,600,800,1000,1200,1400,1600,1800];

tao = nuk/Ek;
expo = -t/tao;

bm = sig0/Em + (sig0/Ek)*(1-exp(expo)) + (sig0/num).*t;

figure(1)
plot(t, bm);
figure(2);
plot(t,fm);

% t0 = 0;
% 
% c1 = 33.3e-5;
% c2 = 1.8391e-10;
% c3 = 5.312e-7;
% c4 = 1.6575e-6;
% str = 5;
% 
% j = c1 + c2*(t-t0) + c3 * (1-exp((t0-t)*(c4/c3)));
% 
% elong = str * j;
% 
% h1 = 31.2;
% matC1 = 4.11e-9;
% matC2 = 1.82;
% matD1 = -11.90;
% matD2 = 1.6776;
% matStr = [10.83 12.48 14.13 16.27 18.32 19.58];
% 
% for k = 1:length(matStr)
%   h2 = log(1+h1) / (matD1 * matStr(k)^(matD2-1));
%   eMat(k,1:length(t)) = (1/(h2*matStr(k))) * log(1 + (1+h1) * exp(h2*matC1*matStr(k)^(matC2+1)*t) - 1);
%   S(k,1:length(t)) = eMat(k)/matStr(k);
%   E(k,1:length(t)) = S(k,1:length(t)).^-1;
% end