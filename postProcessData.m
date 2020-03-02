fileName = 'DOE_Runs';

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

r6M = ((rMax(:,1)-rMax(:,2))./rMax(:,1)) * 100;
r1Y = ((rMax(:,1)-rMax(:,3))./rMax(:,1)) * 100;

rPercLoss = [r6M,r1Y];


