function Q = stiffMat(mstiff, compFunc)
%% -----------------------------------------------------------------------------
% Interpretation of the material stiffness properties vector
% ------------------------------------------------------------------------------
E11 = mstiff(1);
E22 = mstiff(2);
G12 = mstiff(3);
G23 = mstiff(4); % Only used if nu23 is not specified
nu12 = mstiff(5);
nu23 = mstiff(5); % E22/(2*G23)-1;

if ~strcmp(compFunc,'no')
  s = compFunc(mstiff);
  S11 = s(1);
  S22 = s(2);
  S33 = s(3);
  G12 = s(4);
  nu12 = s(5);
  nu23 = s(6);
%   if length(out) == 1
%     S11 = 1/E11;
%     S22 = out;
%     S33 = S22;
%     G12 = 1/G12;
%   elseif length(out) == 2
%     S11 = 1/E11;
%     S22 = out(1);
%     S33 = S22;
%     G12 = out(2); % This was specified for IM&_8552_Tzeng2001
%   elseif length(out) == 3
%     S11 = out(3);
%     S22 = out(1);
%     S33 = out(2);
%     G12 = 1/G12;
%   else
%     error('Unknown compliance function output. Change output, or modify program')
%   end
else
  S11 = 1/E11;
  S22 = 1/E22;
  S33 = S22;
  nu23 = 0.50;
  G12 = 1/G12;
end

%% -----------------------------------------------------------------------------
% Assembling the compliance matrix and inversion to stiffness matrix
% ------------------------------------------------------------------------------
S = [ S11      -nu12*S11 -nu12*S11 0;
     -nu12*S11 S22       -nu23*S33 0;
     -nu12*S11 -nu23*S22 S33       0;
     0         0         0         G12];

Q = inv(S);
% Q = [1 2 3 4; 5 6 7 8; 9 10 11 12; 13 14 15 16];
% Q(1:2:3,1:2:3) = Q(1:2:3,1:2:3)-Q(1:2:3,2:2:4)*...
%     inv(Q(2:2:4,2:2:4))*Q(2:2:4,1:2:3);
