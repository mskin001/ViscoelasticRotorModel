
t = linspace(0,4700,4700);
t0 = zeros(1,length(t));
sigma = 79.25;

c1 = 3.3333e-5;
c2 = 1.8391e-10;
c3 = 5.3121e-7;
c4 = 1.6575e-6;

j = c1 + c2.*(t - t0) + c3*(1 - exp((t0-t)*(c4/c3)));

e = sigma * j;

strain = [e(1), e(100), e(1000), e(2500), e(4700)];
disp(strain)