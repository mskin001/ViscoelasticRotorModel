function [failure, strengthRatio, Fmode, Floc] = failureIndex(Ftype)

global rotor numRims rArr sArr rdiv
% Radial failure is assumed to be rotor damage not totoal failure and is found
% based on the choses failure criterion (Tsai-Wu, MaxR)

if strcmp(Ftype, 'TsaiWu')
%% Tsai-Wu failure criterian
  % Calculates failure index based on Tsai-Wu failure criterion
  for k = 2:numRims
    rvstart = (k-1)*rdiv + 1;
    rvend = k*rdiv;

    LgTens = rotor.stren{k}(1);
    LgComp = rotor.stren{k}(2);
    TranTens = rotor.stren{k}(3);
    TranComp = rotor.stren{k}(4);

    Ftt = 1/(LgTens*LgComp);
    Ft = 1/LgTens - 1/LgComp;
    Frr = 1/(TranTens*TranComp);
    Fr = 1/TranTens - 1/TranComp;
    Ftr = -1/(2*sqrt(LgTens*LgComp*TranTens*TranComp));

    a = Ftt*sArr(1,rvstart:rvend).^2 + 2*Ftr*sArr(1,rvstart:rvend).*sArr(3,rvstart:rvend)...
      + Frr*sArr(3,rvstart:rvend).^2;
    b = Ft*sArr(1,rvstart:rvend) + Fr*sArr(3,rvstart:rvend);

    R1(rvstart:rvend) = (-b + sqrt(b.^2 - 4*a*(-1)))./(2*a);
  end

  strengthRatio = R1.^-1;
  strengthRatio(isinf(strengthRatio)) = [];
  if sum(strengthRatio >= 1)
    failInd = find(strengthRatio >= 1,1);
    radInd = rArr(rdiv+1:end);
    Floc = radInd(failInd);
    
    fpos = [];
    k = 0;
    while isempty(fpos)
      k = k + 1;
      fpos = find((rotor.radii{k}(1) < Floc) && (Floc < rotor.radii{k}(2)),1); % index location of for the outer radius of the failed rim
    end
    hsr = sArr(1,(failInd+rdiv)) / rotor.stren{k}(1); % normalized max hoop stress
    rsr = sArr(3,(failInd+rdiv)) / rotor.stren{k}(3); % normalized max radial stress
    
    if hsr > rsr
      Fmode = 'Burst';
    else
      Fmode = 'Radial';
    end
    failure = 1;
    
  else
    failure = 0;
    Fmode = 'none';
    Floc = 0;
    strengthRatio = 0;
  end

elseif strcmp(Ftype, 'MaxR')
%% Max stress failure criterion
  % Burst failure when the hoop stress exceeds the hoop strength
  for k = 2:numRims
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
  
  % Failure index based on maximum radial Stress
  for k = 2:numRims
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
    Fmode = 'Radial';
    failure = 1;
%     plot(rArr,sArr(3,:));
%     hold on
%     plot(cRadFail, sArr(3,cfi), 'r*')
%     plot(tRadFail, sArr(3,tfi), 'b*')
  else % If no failure occures
    failure = 0;
    Fmode = 'none';
    Floc = 0;
    strengthRatio = 0;
  end

end
