function [k] = failureIndex()

global rotor failure numRims rArr sArr rdiv
% Radial failure is assumed to be rotor damage not totoal failure and is found
% based on the choses failure criterion (Tsai-Wu, MaxR)

if strcmp(failure.crit, 'TsaiWu')
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

    R(rvstart:rvend) = (-b + sqrt(b.^2 - 4*a*(-1)))./(2*a);
  end

  failure.sr = R.^-1;
  failure.sr(isinf(failure.sr)) = [];

  if sum(failure.sr >= 1)
    failInd = find(failure.sr >= 1,1);
    radInd = rArr(rdiv+1:end);
    failure.loc = radInd(failInd);

    fpos = [];
    k = 0;
    while isempty(fpos)
      k = k + 1;
      fpos = find((rotor.radii{k}(1) < failure.loc) && (failure.loc < rotor.radii{k}(2)),1); % index location of for the outer radius of the failed rim
    end
    failure.rim = k;
    hsr = sArr(1,(failInd+rdiv)) / rotor.stren{k}(1); % normalized max hoop stress
    rsr = sArr(3,(failInd+rdiv)) / rotor.stren{k}(3); % normalized max radial stress

    if hsr > rsr
      failure.mode = 'Burst';
    else
      failure.mode = 'Radial';
    end
    failure.occured = true;

  else
    failure.occured = false;
    failure.mode = 'none';
    failure.loc = [];
    failure.sr = [];
  end

elseif strcmp(Ftype, 'MaxR')
%% Max stress failure criterion
  for k = 2:numRims
    rvstart = (k-1)*rdiv + 1;
    rvend = k*rdiv;
    hoopFailureInd(rvstart:rvend) = (sArr(1,rvstart:rvend) / rotor.stren{k}(1)) >= 1; %hoop Failure Index
    cfi(rvstart:rvend) = (sArr(3,rvstart:rvend) / rotor.stren{k}(4)) >= 1;
    tfi(rvstart:rvend) = (sArr(3,rvstart:rvend) / rotor.stren{k}(3)) >= 1;
  end

  if sum(hoopFailureInd) >= 1
    % Burst failure when the hoop stress exceeds the hoop strength
    failure.occured = true;
    failure.mode = 'Burst';

    failure.loc = rArr(hoopFailureInd);
    fpos = [];
    k = 0;
    while isempty(fpos)
      k = k + 1;
      fpos = find((rotor.radii{k}(1) < failure.loc) && (failure.loc < rotor.radii{k}(2)),1); % index location of for the outer radius of the failed rim
    end
    failure.rim = k;

    for k = 2:numRims
      rstart = (k-1)*rdiv + 1;
      rend = k*rdiv;
      failure.sr = sArr(1,rstart:rend) / rotor.stren{k}(1);
    end
  elseif sum(cfi) >= 1 || sum(tfi) >= 1
    % Radial failure when applied stress exceeds tensile strength
    failure.mode = 'Radial';
    failure.occured = true;

    failure.loc = rArr(tfi);
    fpos = [];
    k = 0;
    while isempty(fpos)
      k = k + 1;
      fpos = find((rotor.radii{k}(1) < failure.loc) && (failure.loc < rotor.radii{k}(2)),1); % index location of for the outer radius of the failed rim
    end
    failure.rim = k;

    for k = 2:numRims
      rstart = (k-1)*rdiv + 1;
      rend = k*rdiv;
      failure.sr = sArr(3,rstart:rend) / rotor.stren{k}(3);
    end
  else
    % If failure isn't predicted
    failure.occured = false;
    failure.mode = 'none';
    failure.loc = [];
    failure.sr = [];
    failure.rim = [];
  end
end
