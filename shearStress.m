function [C_shear] = shearStress(alpha, accType, b, w0, tStep, rdiv)

global rim mat tauArr

a = zeros(1:length(rim) - 1);
C = zeros(1:length(rim) - 1);

if strcmp(accType, 'const')
  a(end) = (mat.rho{end} * alpha(tStep)/tStep) / mat.Q{end}(4,4);
else
  a(end) = (mat.rho{end}(rim(end)) * alpha(b*tStep,w0)/(b*tStep)) / mat.Q{end}(4,4);
end

C(end) = (a(end)*(rim(end))^4)/8;

% for k = 1:length(rim) - 2
%   a(end-k) = (mat.rho{end-k} * alpha(tStep,w0)/tStep) / mat.Q{end}(4,4);
%   gRatio = mat.Q{end-k+1}(4,4)/mat.Q{end-k}(4,4);
%   C(end-k) = (a(end-k) - gRatio * a(end-k+1)) * rim(k+1)^4/8 + gRatio*C(end-k+1);
% end
C_shear = C;

for k = 1:length(rim) - 1 
    dr = linspace(rim(k),rim(k+1),rdiv);
    rStart = (k-1)*rdiv + 1;
    rEnd = k*rdiv;
    tauArr(1,rStart:rEnd) = (-mat.Q{1,k}(4,4) * ((a(k)*dr.^2)./4 - (2*C(k))./dr.^2));
end

