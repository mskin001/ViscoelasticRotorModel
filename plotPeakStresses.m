resFile = 'constResults.mat';
load(resFile)
load('MaterialProperties\CFRP_Aparicio2011.mat');

figure()
hold on
maxStr = zeros(length(results.time),3);
for k = 1:length(results.time)
  hoopStr = results.sArr{k}(1,:,1);
  radStr = results.sArr{k}(3,:,1);
  shear = results.tauArr{k};
  maxStr(k,:) = [max(hoopStr), max(radStr), max(shear)];
  
end

nStr = [maxStr(:,1)/stren(1), maxStr(:,2)/stren(3), maxStr(:,3)/stren(end)];

figure(1)
hold on
plot(results.vel,maxStr(:,1)*10^-6, '-*', 'Color', [0 0.4470 0.7410], 'MarkerIndices', 1:10:length(results.time), 'Linewidth', 1.5)
plot(results.vel,maxStr(:,2)*10^-6, '--d', 'Color', [0.6350 0.0780 0.1840], 'MarkerIndices', 1:10:length(results.time), 'Linewidth', 1.5)
% plot(results.time,maxStr(:,3)*10^-6, 'Linewidth', 1.5)

ylabel('Maximum stress [MPa]')
xlabel('Angular velocity [rpm]')
legend('Circ. Stress', 'Radial Stress', 'Location', 'northwest')
set(gca, 'Fontsize', 12)

figure(2)
hold on
plot(results.vel,nStr(:,1), '-*', 'Color', [0 0.4470 0.7410], 'MarkerIndices', 1:10:length(results.time), 'Linewidth', 1.5)
plot(results.vel,nStr(:,2), '--d', 'Color', [0.6350 0.0780 0.1840], 'MarkerIndices', 1:10:length(results.time), 'Linewidth', 1.5)
plot(results.vel,results.peakstr, ':v', 'Color', [0.4660 0.6740 0.1880], 'MarkerIndices', 1:10:length(results.time), 'Linewidth', 1.5)
% plot(results.time,nStr(:,3), 'Linewidth', 1.5)

ylabel('Normalzied stress')
xlabel('Angular velocity [rpm]')
legend('Circ. Stress', 'Radial Stress', 'Tsai-Wu SR', 'Location', 'northwest')
set(gca, 'Fontsize', 12)




% rArr = linspace(0.08, 0.2,30);
% figure(1)
% hold on
% plot(rArr,results.sArr{1}(3,:,1)*10, 'b-', 'Linewidth', 2)
% plot(rArr,results.sArr{1}(1,:,1), 'r-.', 'Linewidth', 2)
% A2011 = csvread('Apa2011_10xradial.csv');
% 
% radii = A2011(:,1);
% str = A2011(:,2);
% figure(1)
% hold on
% plot(radii,str,'k*')
% 
% 
% A2011 = csvread('Apa2011_hoop.csv');
% 
% radii = A2011(:,1);
% str = A2011(:,2);
% figure(1)
% hold on
% plot(radii,str,'kd')
% 
% xlabel('Radius [m]')
% ylabel('Stress [MPa]')
% legend('Model 10*\sigma_r', 'Model \sigma_\theta', 'Aparicio 10*\sigma_r', 'Aparicio \sigma_\theta')
% set(gca, 'Fontsize', 12)