Et = 38.6e9;
Ez = 8.27e9;
Er = Ez;

Gtz = 4.14e9;
Gtr = Gtz;
Grz = Gtz;

nutz = 0.26;
nutr = nutz;
nurz = nutz;

Qha = [1/Et -nutr/Er -nutz/Ez;
       -nutr/Er 1/Er -nurz/Ez;
       -nutz/Ez -nurz/Ez 1/Ez]^-1
     
     
load('MaterialProperties/Glass_Epoxy_Ha1999.mat')
E11 = mstiff(1);
E22 = mstiff(2);
G12 = mstiff(3);
G23 = mstiff(4); % Only used if nu23 is not specified
nu12 = mstiff(5);
nu23 = nu12; % E22/(2*G23)-1;

S11 = 1/E11;
S22 = 1/E22;
S33 = S22;
G12 = 1/G12;

S = [ S11    -nu12*S22 -nu12*S33 0;
     -nu12*S22 S22     -nu23*S33 0;
     -nu12*S33 -nu23*S22  S33 0;
     0 0 0 G12];
Q = inv(S);
% Q(1:2:3,1:2:3) = Q(1:2:3,1:2:3)-Q(1:2:3,2:2:4)*...
%     inv(Q(2:2:4,2:2:4))*Q(2:2:4,1:2:3)