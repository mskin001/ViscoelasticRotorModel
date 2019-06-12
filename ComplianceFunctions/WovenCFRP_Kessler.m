function [varargout] = WovenCFRP_Kessler()
% This fuction produces the time dependent compliance of the woven carbon
% fiber epoxy composite studied by Goertzen and Kessler in "Creep behavior
% of carbon fiber/epoxy composites." The compliance master curve was taken
% from figure 11a useing a pdf to graph converter and curve fitting the
% discrete data using matlab's built in curve fitting toolbox. The data was
% fit uings a one term exponential function with an r^2 value greater than
% 0.98.
%
% This function uses the total time to predict the compliance. This is
% taken as the S22, compliance as of writing this, but note the material is
% woven so this is not entirely correct.

global t

a = 2.181e-10;
b = 0.1122;

varargout{1} = a .* exp(b .* log10(t));
end
