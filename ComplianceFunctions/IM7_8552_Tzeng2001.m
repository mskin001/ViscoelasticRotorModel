function [s] = IM7_8552_Tzeng2001(mstiff)
% This function reproduces the compliance curve for IM7/8552 discussed by
% Tzeng 2001 in "Viscoelastic Analysis of Composite Cylinders Subjected to
% Rotation." The compliance was as a best fit equation in eq.19 of this
% paper. 
%
% The units given in the paper are inches and pounds. The compliance is
% in^2/lb

global t

tr = mstiff(2);
sh = mstiff(3);

s(1) = 1/mstiff(1);
s(2) = 1/tr * (t/3600)^0.03;
s(3) = s(2);
s(4) = 1/sh * (t/3600)^0.03;
s(5) = mstiff(5);
s(6) = 0.36;