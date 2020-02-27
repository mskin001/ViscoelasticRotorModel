function [s] = Militky15(mstiff)

global t

Em = mstiff(1);
Ek = mstiff(2);
eta_m = mstiff(3);
eta_k = mstiff(4);

tao = eta_k/Ek;
expo = -t/(tao);

s(1:4) = 1/Em + (1/Ek)*(1-exp(expo)) + t/eta_m;
s(5:6) = mstiff(end);