resFile = 'constResults.mat';
load(resFile)
load('MaterialProperties\CFRP_Aparicio2011.mat');

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
plot(results.time,maxStr(:,1)*10^-6, 'Linewidth', 1.5)
plot(results.time,maxStr(:,2)*10^-6, 'Linewidth', 1.5)
plot(results.time,maxStr(:,3)*10^-6, 'Linewidth', 1.5)

ylabel('Maximum stress [MPa]')
xlabel('Time [s]')
legend('Circ. Stress', 'Radial Stress', 'Shear Stress', 'Location', 'northwest')
set(gca, 'Fontsize', 12)

figure(2)
hold on
plot(results.time,nStr(:,1), 'Linewidth', 1.5)
plot(results.time,nStr(:,2), 'Linewidth', 1.5)
plot(results.time,nStr(:,3), 'Linewidth', 1.5)

ylabel('Normalzied stress [MPa]')
xlabel('Time [s]')
legend('Circ. Stress', 'Radial Stress', 'Shear Stress', 'Location', 'northwest')
set(gca, 'Fontsize', 12)
