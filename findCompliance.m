fileName = 'IM7 creep strain.csv';
stress = 41.37e6; % Pa

rawData = csvread(fileName);
xData = rawData(:,1); % time
yData = rawData(:,2); % strain

ind = find(xData,1,'first');

xData(1:ind-1) = [];
yData(1:ind-1) = [];
comp = yData / stress;

a = 1.395e-10;
b = 2.912e-7;
c = -7.472e-12;
d = -3.708e-5;
t = linspace(1,250000,250000);

s22 = a.*exp(b.*t) + c.*exp(d.*t);

figure(1)
hold on
plot(xData,comp, 'bo')
plot(t,s22,'r-')

% figure(2)
% plot(log10(xData),comp)