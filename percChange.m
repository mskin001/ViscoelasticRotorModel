function [] = percChange(material)
global sArr

radStr_t1 = abs(sArr(3,:,1));
radStr_tend = abs(sArr(3,:,end));

[max_t1, ind] = max(radStr_t1);
max_tend = radStr_tend(ind);

radStr_dif = max_t1 - max_tend;
radStr_perc = (radStr_dif/max_t1) * 100;
radStrResult = [radStr_dif, radStr_perc];

hoopStr_t1 = abs(sArr(1,:,1));
hoopStr_tend = abs(sArr(1,:,end));
[max_t1, ind] = max(hoopStr_t1);
max_tend = hoopStr_tend(ind);

hoopDif = hoopStr_t1 - max_tend;
hoopPerc = (hoopDif/max_t1) * 100;
