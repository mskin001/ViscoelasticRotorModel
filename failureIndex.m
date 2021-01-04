function [strRatio] = failureIndex(rdiv)

global results rArr rim mat vari

for b = 1:vari+1
  for k = 1:length(rim) - 1
    rStart = (k-1)*rdiv + 1;
    rEnd = k*rdiv;

    LgTens = mat.stren{k}(1); % longitudinal tensile strength
    LgComp = mat.stren{k}(2); % longitudinal compressive strength
    TranTens = mat.stren{k}(3); % transverse tensile strength
    TranComp = mat.stren{k}(4); % transverse compressive strength
    sh = mat.stren{k}(5); % longitudinal shear strength

    Ftt = 1/(LgTens*LgComp); % F11
    Ft = 1/LgTens - 1/LgComp; % F1
    Frr = 1/(TranTens*TranComp); % F22
    Fr = 1/TranTens - 1/TranComp; % F2
    Ftr = -1/(2*sqrt(LgTens*LgComp*TranTens*TranComp)); % F12
    Fs = 1/sh^2; % F66

    sigt = results.sArr{b}(1,rStart:rEnd); % sig1
    sigr = results.sArr{b}(3,rStart:rEnd); % sig2 = sig3
    tau = results.tauArr{b}(1,rStart:rEnd); % tau12

    % F(b, rStart:rEnd) = Ftt*sigt.^2 + 2*Ftr*sigt.*sigr + Frr*sigr.^2 + Fs*tau.^2 + Ft*sigt + Fr*sigr;
    A = Ftt*sigt.^2 + 2*Ftr*sigt.*sigr + Frr*sigr.^2 + Fs*tau.^2;
    B = Ft*sigt + Fr*sigr;
    C = -1;
    R(b,rStart:rEnd) = (-B + sqrt(B.^2 - 4*A*C)) ./ (2*A);
    strRatio(b,rStart:rEnd) = R(b,rStart:rEnd).^-1;

  end
  [peakStr(b), ind] = max(strRatio(b,:));
  peakLoc(b) = rArr(ind);
end

results.peakstr = peakStr;
results.peakloc = peakLoc;
