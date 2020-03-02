fileName = 'DOE_Runs';

addpath('DOE_Runs')
files = dir(fileName);

for k = 3:length(files)
  load(files(k).name);
  hoop = sArr(1,:,:);
  rad = sArr(3,:,:);
  
%   for i = 1:3
%     hMax(i,(k-2)) = max(abs(hoop(:,:,i)));
%     rMax(i,(k-2)) = max(abs(rad(:,:,i)));
%   end

  figure(1)
  hold on
  plot(hoop(:,:,1))
  plot(hoop(:,:,2))
  plot(hoop(:,:,3))
  
  figure(2)
  hold on
  plot(rad(:,:,1))
  plot(rad(:,:,2))
  plot(rad(:,:,3))

end

hDiff = diff(hMax);
hPerc = (-100 * (hDiff./hMax(1:2,:)))';
h6M = vec2mat(hPerc(:,1),3);
h1Y = vec2mat(hPerc(:,2),3);

rDiff = diff(rMax);
rPerc = (-100 * (rDiff./rMax(1:2,:)))';
r6M = vec2mat(rPerc(:,1),3);
r1Y = vec2mat(rPerc(:,2),3);

