function [s] = CFRP_BM_Almeida2017(mstiff)

global t

Em = mstiff(1);
Ek = mstiff(2);
num = mstiff(3);
nuk = mstiff(4);

tao = nuk/Ek;
expo = -t/tao;

s(1) = 1/Em + (1/Ek)*(1-exp(expo)) + (1/num)*t;

Em = mstiff(5);
Ek = mstiff(6);
num = mstiff(7);
nuk = mstiff(8);

tao = nuk/Ek;
expo = -t/tao;

s(2) = 1/Em + (1/Ek)*(1-exp(expo)) + (1/num)*t;

nu12 = mstiff(9);
nu23 = nu12;

s(3:4) = s(2);
s(5:6) = mstiff(9);
