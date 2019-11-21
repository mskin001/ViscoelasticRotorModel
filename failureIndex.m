function [Fmode, failureIndex, Floc] = failureIndex(Ftype)

global mat rArr sArr


%% -----------------------------------------------------------------------------
%  Determine if and burst failure occurs
%  -----------------------------------------------------------------------------
%  Burst failure occurs when the hoop stress exceeds the hoop strength

hoopFailureInd = find((sArr(1,:,1) - mat.stren{1}(1)) >= 1); %hoop Failure Index
if hoopFailure >= 1
  Fmode = 'Burst';
  Floc = rArr(hoopFailureInd);
  return
end


%% -----------------------------------------------------------------------------
%  Determine if radial failure occurs
%  -----------------------------------------------------------------------------
% Radial failure is assumed to be rotor damage not totoal failure and is found
% based on the choses failure criterion (Tsai-Wu, MaxR)

Fmode = 'Radial';

if strcmp(Ftype, 'TsaiWu')
  % Calculates failure index based on Tsai-Wu failure criterion
  LgTens = mat.stren{1}(1);
  LgComp = mat.stren{1}(2);
  TranTens = mat.stren{1}(3);
  TranComp = mat.stren{1}(4);

  Ftt = 1/(LgTens*LgComp);
  Ft = 1/LgTens - 1/LgComp;
  Frr = 1/(TranTens*TranComp);
  Fr = 1/TranTens - 1/TranComp;
  Ftr = -1/(2*sqrt(LgTens*LgComp*TranTens*TranComp));

  a = Ftt*sArr(1,:,1).^2 + 2*Ftr*sArr(1,:,1).*sArr(3,:,1) + Frr*sArr(3,:,1).^2;
  b = Ft*sArr(1,:,1) + Fr*sArr(3,:,1);

  R1 = (-b + sqrt(b.^2 - 4*a.*c))/(2*a);
  R2 = (-b - sqrt(b.^2 - 4*a.*c))/(2*a);

  if all(R1 >= 0)
    failureIndex = R1;
    disp('failureIndex is R1');
  elseif all(R2 >= 0)
    failureIndex = R2;
    disp('failureIndex is R2');
  end

  Floc = rArr(failureIndex);

elseif strcmp(Ftype, 'MaxR')
  % Failure index based on maximum radial Stress
  Cfailure = sArr(3,:,1) - mat.stren{1}(4);
  Tfailure = sArr(3,:,1) - mat.stren{1}(3);

  CFindex = find(cFailure >= 1);
  TFindex = find(Tfailure >=1);

  failureIndex = [CFindex, TfIndex];
  Floc = rArr(failureIndex);
end
