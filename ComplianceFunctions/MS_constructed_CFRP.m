function [s] = MS_constructed_CFRP(mstiff)

global t

Em = mstiff(2);
Ek = mstiff(3);
nu_m = mstiff(4);
nu_k = mstiff(5);

tau = nu_k/Ek;
expo = -t/tau;

s(1) = 1/mstiff(1); % fiber direction
s(2:4) = 1/Em + (1/Ek)*(1-exp(expo)) + t/nu_m; % transverse
% s(4) = (1 + mstiff(6)) / mstiff(1); % shear
s(5:6) = mstiff(6); %poisons ratio
end
