syms ri ro fi1 fi2 fi9 fi10 C1 C2 k

Is = [-1 0; 0 1]
G = [ri^(k+1) ri^(-k+1); ro^(k+1) ro^(-k+1)]
io = [fi1 0; 0 fi2]
iz = [fi9 0; 0 fi10]
C = [C1; C2]
K = [1/(k+1) 0; 0 -1/(k+1)];

% haU = G*io*C
test = Is*G*K*iz*C