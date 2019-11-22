function [failure, strengthRatio, Fmode, Floc] = failureIndex(Ftype)

global rotor numRims rArr sArr rdiv


%% -----------------------------------------------------------------------------
%  Determine if and burst failure occurs
%  -----------------------------------------------------------------------------
%  Burst failure occurs when the hoop stress exceeds the hoop strength
for k = 1:numRims
  rvstart = (k-1)*rdiv + 1;
  rvend = k*rdiv;
  hoopFailureInd(rvstart:rvend) = (sArr(1,rvstart:rvend) / rotor.stren{k}(1)) >= 1; %hoop Failure Index
end

if sum(hoopFailureInd) >= 1
  failure = true;
  Fmode = 'Burst';
  Floc = rArr(hoopFailureInd);
  strengthRatio = sArr(1,rvstart:rvend) / rotor.stren{k}(1);
  fprintf('Burst Failure\n')
  plot(rArr,sArr(1,:));
  hold on
  plot(Floc,sArr(1,hoopFailureInd), 'r*')
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
  for k = 1:length(rim) - 1
    rvstart = (k-1)*rdiv + 1;
    rvend = k*rdiv;

    LgTens = mat.stren{k}(1);
    LgComp = mat.stren{k}(2);
    TranTens = mat.stren{k}(3);
    TranComp = mat.stren{k}(4);

    Ftt = 1/(LgTens*LgComp);
    Ft = 1/LgTens - 1/LgComp;
    Frr = 1/(TranTens*TranComp);
    Fr = 1/TranTens - 1/TranComp;
    Ftr = -1/(2*sqrt(LgTens*LgComp*TranTens*TranComp));

    a = Ftt*sArr(1,rvstart:rvend,1).^2 + 2*Ftr*sArr(1,rvstart:rvend,1).*sArr(3,rvstart:rvend,1) + Frr*sArr(3,rvstart:rvend,1).^2;
    b = Ft*sArr(1,rvstart:rvend,1) + Fr*sArr(3,rvstart:rvend,1);

    R1(rvstart:rvend) = (-b + sqrt(b.^2 - 4*a.*c))/(2*a);
    R2(rvstart:rvend) = (-b - sqrt(b.^2 - 4*a.*c))/(2*a);
  end

  if all(R1 >= 0)
    strengthRatio = R1;
    disp('failureIndex is R1');
  elseif all(R2 >= 0)
    strengthRatio = R2;
    disp('failureIndex is R2');
  end

  Floc = rArr(strengthRatio);

elseif strcmp(Ftype, 'MaxR')
  % Failure index based on maximum radial Stress
  for k = 1:numRims
    rvstart = (k-1)*rdiv + 1;
    rvend = k*rdiv;
    cfi(rvstart:rvend) = (sArr(3,rvstart:rvend) / rotor.stren{k}(4)) >= 1;
    tfi(rvstart:rvend) = (sArr(3,rvstart:rvend) / rotor.stren{k}(3)) >= 1;
  end

  if sum(cfi) >= 1 || sum(tfi) >= 1 % If there is a tensile or compressive radial failure
    cRadFail = rArr(cfi);
    tRadFail = rArr(tfi);
    strengthRatio = 1; % Still need to work on this. Compressive or tensile failure? Both? Under construction
    Floc = [cRadFail, tRadFail];
    failure = true;
    plot(rArr,sArr(3,:));
    hold on
    plot(cRadFail, sArr(3,cfi), 'r*')
    plot(tRadFail, sArr(3,tfi), 'b*')
  else % If no failure occures
    failure = false;
    Fmode = 'none';
    Floc = 0;
    strengthRatio = 0;
  end

end
