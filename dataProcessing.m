% clear variables
fileName = 'DOE_Runs';

addpath('DOE_Runs')
files = dir(fileName);

for k = 3:length(files)
  load(files(k).name);
  hoop = result.sArr(1,:,:);
  rad = result.sArr(3,:,:);
  
  for i = 1:3
    hMax(i,(k-2)) = max(abs(hoop(:,:,i)));
    rMax(i,(k-2)) = max(abs(rad(:,:,i)));
  end
% 
%   figure()
%   hold on
%   plot(hoop(:,:,1))
%   plot(hoop(:,:,2))
%   plot(hoop(:,:,3))
%   
%   figure()
%   hold on
%   plot(rad(:,:,1))
%   plot(rad(:,:,2))
%   plot(rad(:,:,3))

end

hDiff = diff(hMax);
hPerc = (-100 * (hDiff./hMax(1:2,:)))';
h1Y = vec2mat(hPerc(:,1),3);
h5Y = vec2mat(hPerc(:,2),3);

rDiff = diff(rMax);
rPerc = (-100 * (rDiff./rMax(1:2,:)))';
r1Y = vec2mat(rPerc(:,1),3);
r5Y = vec2mat(rPerc(:,2),3);

mr_nuM(1,:) = -1*sum(h1Y)/3; %C1Y
mr_nuM(2,:) = -1*sum(h5Y)/3; %C5Y
mr_nuM(3,:) = sum(r1Y)/3; %R1Y
mr_nuM(4,:) = sum(r5Y)/3; %R5Y

mr_nuK(1,:) = -1*sum(h1Y')/3;
mr_nuK(2,:) = -1*sum(h5Y')/3;
mr_nuK(3,:) = sum(r1Y')/3;
mr_nuK(4,:) = sum(r5Y')/3;
xVal = [1 2 3];

hold on
figure(1)
plot(xVal, mr_nuM(1,:), 'bo--', 'LineWidth', 2)
plot(xVal, mr_nuM(2,:), 'o-', 'Color', [0 0.4470 0.7410], 'LineWidth', 2)
plot(xVal, mr_nuM(3,:), 'r:^', 'LineWidth', 2)
plot(xVal, mr_nuM(4,:), '-.', 'Color', [0.6350 0.0780 0.1840], 'LineWidth', 2)
legend('Cir. Str 1Y', 'Cir. Str 5Y', 'Rad. Str 1Y', 'Rad. Str 5Y', 'Location', 'southoutside', 'NumColumns', 4)
xticks([1 2 3])
xlabel('\it \eta_M \rm level')
ylabel('Mean Response')
set(gca, 'FontSize', 12)
figure(2)
hold on
plot(xVal, mr_nuK(1,:), 'bo--', 'LineWidth', 2)
plot(xVal, mr_nuK(2,:), 'o-', 'Color', [0 0.4470 0.7410], 'LineWidth', 2)
plot(xVal, mr_nuK(3,:), 'r:^', 'LineWidth', 2)
plot(xVal, mr_nuK(4,:), '-.', 'Color', [0.6350 0.0780 0.1840], 'LineWidth', 2)
legend('Cir. Str 1Y', 'Cir. Str 5Y', 'Rad. Str 1Y', 'Rad. Str 5Y', 'Location', 'southoutside',...
  'NumColumns',4)
xticks([1 2 3])
xlabel('\it \eta_K \rm level')
ylabel('Mean Response')
set(gca, 'FontSize', 12)

  
  
  
