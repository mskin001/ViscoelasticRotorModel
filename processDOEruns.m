clear variables; clc

addpath('DOE_Runs')

for k = 1:9
  fileName = ['Run', num2str(k), '.mat'];
  load(fileName);
  hoopStr = result.sArr(1,:,:);
  radStr = result.sArr(3,:,:);

  for b = 1:3
    maxVal.hoop(k,b) = max(hoopStr(1,:,b));
    maxVal.rad(k,b) = max(abs(radStr(1,:,b)));
    if b == 1
      del.Hoop(k,b) = 100 * ((maxVal.hoop(k,b) - maxVal.hoop(k,1)) / maxVal.hoop(k,1));
      del.Rad(k,b) = 100 * ((maxVal.rad(k,1) - maxVal.rad(k,b)) / maxVal.rad(k,1));
    else
      del.Hoop(k,b) = 100 * ((maxVal.hoop(k,b) - maxVal.hoop(k,b-1)) / maxVal.hoop(k,b-1));
      del.Rad(k,b) = 100 * ((maxVal.rad(k,b-1) - maxVal.rad(k,b)) / maxVal.rad(k,b-1));
    end
    
  end
end


