function [tau] = shearStress(alpha, rdiv)

global rim mat

a = zeros(1:length(rim) - 1);
C = zeros(1:length(rim) - 1);

a(end) = (mat.rho{end} * alpha) / mat.Q{end}(4,4);
C(end) = (a(end)*rim(end)^4)/8;

b = 1;

for k = 1:length(rim) - 2
  a(end-k) = (mat.rho{end-k} * alpha) / mat.Q{end}(4,4);
  gRatio = mat.Q{end-k+1}(4,4)/mat.Q{end-k}(4,4);
  C(end-k) = (a(end-k) - gRatio * a(end-k+1)) * rim(k+1)^4/8 + gRatio*C(end-k+1);
end

for k = 1:length(rim) - 1 
    dr = linspace(rim(k),rim(k+1),rdiv);
    rStart = (k-1)*rdiv + 1;
    rEnd = k*rdiv;
    tau(1,rStart:rEnd) = mat.Q{b,k}(4,4) * (-(a(k)*dr.^2)./4 + (2*C(k))./dr.^2);
end

disp(tau)