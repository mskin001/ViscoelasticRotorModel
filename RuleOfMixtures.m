% Ef = 240e9;
% Em = 2.419e9;
% f = 0.72;
% 
% E11 = f*Ef + (1-f)*Em;
% E22 = (f/Ef + (1-f)/Em)^-1;

t = linspace(1,5*8760,5*8760);

mstiff = [1.7348e11, 8.4e9, 1.173e9, 8.063e13, 1.715e12, 0.3];

Em = mstiff(2);
Ek = mstiff(3);
nu_m = mstiff(4);
nu_k = mstiff(5);

tau = nu_k/Ek;
expo = -t/tau;

% s(1) = 1/mstiff(1); % fiber direction
s = 1/Em + (1/Ek)*(1-exp(expo)) + t/nu_m; % transverse
% s(4) = (1 + mstiff(6)) / mstiff(1); % shear
% s(5:6) = mstiff(6); %poisons ratio

plot(t,s*10^9, 'LineWidth', 1.5)
xlabel('time [hr]')
ylabel('Compliance [1/GPa]')
set(gca,'FontSize',12)