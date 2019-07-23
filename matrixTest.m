syms ri ro fi1 fi2 fi9 fi10 C1 C2 k


G = [-ri^(k+1) -ri^(-k+1); ro^(k+1) ro^(-k+1)]
iz = [fi9/(k+1) 0; 0 (fi10/(-k+1))]
C = [C1; C2]


% haU = G*io*C
test0 = iz*G
test1 = test0*C

test = G*iz*C