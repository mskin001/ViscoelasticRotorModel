function [s] = Almeida2018(mstiff)

global t

Em = mstiff(1);
Ek = mstiff(2);
eta_m = mstiff(3);
eta_k = mstiff(4);

tao = eta_k/Ek;
expo = -t/(tao);

s(1) = 1/Em + (1/Ek)*(1-exp(expo)) + t/eta_m;

Em = mstiff(5);
Ek = mstiff(6);
eta_m = mstiff(7);
eta_k = mstiff(8);

tao = eta_k/Ek;
expo = -t/tao;

s(2) = 1/Em + (1/Ek)*(1-exp(expo)) + t/eta_m;
s(3:4) = s(2);
s(5:6) = mstiff(9);
