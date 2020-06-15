t = linspace(0,10,100);
b = 1;
tau = -0.5;

expo = b*exp(tau*t);
hold on
plot(t,expo,'LineWidth',2)

% t = linspace(0,pi/4,100);
% b = -2; % amplitude
% x = 4; % period
% c = (3*pi)/2; % horizontal phase shift
% d = b; % vertical phase shift
% sinu = b*sin(x*t + c)+d;


plot(t,expo,'LineWidth',2)

% figure(2)
% plot(t,sinu,'LineWidth',2)