fileName = 'DOE_Runs';
numLevels = 3;
factorLevels = [1 1; 2 1; 3 1; 1 2; 2 2; 3 2; 1 3; 2 3; 3 3];
levels = [1,2,3];

addpath(fileName);
files = dir(fileName);

% b = 1;
% l = 2;
for k = 3:length(files)

    load(files(k).name);
    rad = sArr(3,:,:);
    hoop = sArr(1,:,:);

    for i = 1:3
        [rMax((k-2),i),indR] = max(abs(rad(:,:,i)));
        [hMax((k-2),i),indH] = max((hoop(:,:,i)));

%         figure(b), hold on
%         plot(rad(1,:,i),'LineWidth',1.5);
%         plot(indR,-1*rMax((k-2),i),'bo');
%
%         figure(l), hold on
%         plot(hoop(1,:,i),'LineWidth',1.5);
%         plot(indR,hMax((k-2),i),'k^');
    end
%
%     b = b + 2;
%
%     l = l + 2;
end

r6M = ((rMax(:,1)-rMax(:,2)) ./ rMax(:,1)) * 100;
r1Y = ((rMax(:,1)-rMax(:,3)) ./ rMax(:,1)) * 100;
rpl = [r6M, r1Y];

h6M = ((hMax(:,1) - hMax(:,2)) ./ hMax(:,1)) * 100;
h1Y = ((hMax(:,1) - hMax(:,3)) ./ hMax(:,1)) * 100;
hpl = [h6M, h1Y];


for k = 1:3
  mLevelInd = find(factorLevels(:,1) == k); %find ind for eta_m factor levels
  kLevelInd = find(factorLevels(:,2) == k); %find ind for eta_k factor levels

  mR_resp = rpl(mLevelInd,:);
  mH_resp = hpl(mLevelInd,:);

  kR_resp = rpl(kLevelInd,:);
  kH_resp = hpl(kLevelInd,:);

  mR_eff(k,1:2) = sum(mR_resp) ./ numLevels;
  mH_eff(k,1:2) = sum(mH_resp) ./ numLevels;
  kR_eff(k,1:2) = sum(kR_resp) ./ numLevels;
  kH_eff(k,1:2) = sum(kH_resp) ./ numLevels;
end


figure(1); hold on
plot(levels, mR_eff(:,1), 'b-o', 'LineWidth', 2);
plot(levels, mR_eff(:,2), '--o', 'Color', [0 0.4470 0.7410], 'LineWidth', 2);
plot(levels, mH_eff(:,1), 'r-.^', 'LineWidth', 2);
plot(levels, mH_eff(:,2), ':^', 'Color', [0.6350 0.0780 0.1840], 'LineWidth', 2);
xticks([1 2 3])
xlabel('\eta_m level')
ylabel('Mean Response')
legend('Rad Str 6M', 'Rad Str 1Y', 'Hoop Str 6M', 'Hoop Str 1Y')
set(gca, 'FontSize', 12)

figure(2); hold on
plot(levels, kR_eff(:,1), 'b-o', 'LineWidth', 2);
plot(levels, kR_eff(:,2), '--o', 'Color', [0 0.4470 0.7410], 'LineWidth', 2);
plot(levels, kH_eff(:,1), 'r-.^', 'LineWidth', 2);
plot(levels, kH_eff(:,2), ':^', 'Color', [0.6350 0.0780 0.1840], 'LineWidth', 2);
xticks([1 2 3])
xlabel('\eta_k level')
ylabel('Mean Response')
legend('Rad Str 6M', 'Rad Str 1Y', 'Hoop Str 6M', 'Hoop Str 1Y', 'Location', 'southeast')
set(gca, 'FontSize', 12)
