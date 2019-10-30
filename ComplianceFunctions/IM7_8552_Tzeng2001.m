function [out] = IM7_8552_Tzeng2001(mstiff)
% This function reproduces the compliance curve for IM7/8552 discussed by
% Tzeng 2001 in "Viscoelastic Analysis of Composite Cylinders Subjected to
% Rotation." The compliance was as a best fit equation in eq.19 of this
% paper. 
%
% The units given in the paper are inches and pounds. The compliance is
% in^2/lb

global t

tr = mstiff(2);
s = mstiff(3);

out(1) = 1/tr * (t)^0.03;
out(2) = 1/s * (t)^0.03;