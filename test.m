t = linspace(0,87600,87600);
A = 5.43e-2;
n = 0.0117;
e0 = 0;
fm = e0 + A*t.^n;

Em = 1.69e14;
num = 1.03e13;
Ek = 1636.76e6;
nuk = 1879.26e6;
sig0 = 1;

% sig0 = 1;
% Em = 2.096e9;
% num = 1.67e13;
% Ek = 17.76e9;
% nuk = 1.32e12;
% t = [0,10,20,40,60,80,100,200,300];

tau = nuk/Ek;
expo = -t/tau;

fiber = sig0*(1/Em + (1/Ek)*(1-exp(expo)) + t/num);

% sig0 = 1;
Em = 6.75e13;
num = 3.60e12;
Ek = 171.76e6;
nuk = 170.73e6;

tau = nuk/Ek;
expo = -t/tau;

trans = sig0*(1/Em + (1/Ek)*(1-exp(expo)) + t/num);

time = log10(t);

figure(1);
hold on
plot(time,trans*10^9,'b-', 'LineWidth', 1.5)
plot(time,fiber*10^9, 'r--', 'LineWidth',1.5)
xlabel('Log(t) [hr]')
ylabel('Compliance [1/GPa]')
set(gca, 'FontSize',12)
legend('Transverse - [90]_{12}', 'Fiber - [0]_{12}', 'Location', 'northwest')

% figure(2)
% plot(t,bm)