rawData = csvread('IM7_8552_Tzeng2001.csv');
t = rawData(2:end,1);
compliance = rawData(2:end,2);


s = 7.53284e-7;

y = log10(s * log10(t).^0.03);

% y = a * exp(b*t) + c * exp(d*t);
% plot(t,compliance, 'bo');
hold on
plot(t,y, 'r-')