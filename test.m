rpm = 50000;
rpmMax = 98000;
tStep = 0.2;
tMax = 20;
T = tMax / log(rpmMax/rpm);

w0 = rpm * pi / 30;

for k = 1:100
  time(k) = k*tStep;
  expo(k) = w0*exp((k*tStep)/T);
end
hold on
plot(time,expo,'LineWidth',2)

% t = linspace(0,20,101);
% t = 0.2;
% b = 2356.2; % amplitude
% x = (2*pi/40); % period
% c = (3*pi)/2; % horizontal phase shift
% d = b; % vertical phase shift
% 
% for k = 1:100
%   time(k) = t*k;
%   sinu(k) = b*sin((x*k*t) + c)+d;
%   vel(k) = 50000*pi/30 + sinu(k);
% end
% % figure(2)
% plot(time,sinu,'LineWidth',2)
% hold on
% plot(time,vel,'LineWidth',2)