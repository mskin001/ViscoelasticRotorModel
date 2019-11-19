function [C] = discretizeStressStrain(rdiv, delta)
% This function calculates a discrete vector of stress and strain in each rim of
% the rotor. It uses the boundary conditions for each rim caulculated in the
% function boundaryConditions to find the C constants. These are then used to
% find a discrete vector of displacements in that rim. The displacements are
% converted to strain then stress for that rim. The displacements, stress, and
% strain are stored in their respective arrays.
%
% This function was designed for use with the creep strain model designed by
% Miles Skinner
%
% This code was created based on the model described by Ha et al. in "Optimmum
% Design of Thick-Walled Composite Rings for an Energy Storage System," and
% was informed by the optimization code created by Krack (04/2009).
%
% Programmer: Miles Skinner, Advanced Composite Materials Engineering,
%             University of Alberta
% Origional Language: MATLAB 2018b

%% -----------------------------------------------------------------------------
% Preallocate variables
%-------------------------------------------------------------------------------
global mat rim U w rArr uArr sArr eArr vari

%% -----------------------------------------------------------------------------
% Calculate stress, strain, displacement for each rim
% ------------------------------------------------------------------------------
prog = waitbar(0,'Descretize Stress/Strain', 'CreateCancelBtn',...
  'setappdata(gcbf,''Canceling'',1)');
setappdata(prog,'Canceling',0);

for k = 1:length(rim)-1
  Q11 = mat.Q{b,k}(1,1);
  Q12 = mat.Q{b,k}(1,2);
  Q13 = mat.Q{b,k}(1,3);
  Q23 = mat.Q{b,k}(2,3);
  Q33 = mat.Q{b,k}(3,3);
  kappa = sqrt(Q11/Q33); % intermediate variable of stiffness ratio

  % intermediate variables for calculating the stiffness matrix
  fi0 = 1/(Q33 * (9 - kappa^2));
  fi1 = 1/(Q13 + (kappa * Q33));
  fi2 = 1/(Q13 - (kappa * Q33));

  % Calculate absolute displacement at the inner and outer surface of each rim
  u = [U(b,k); U(b,k+1) - delta(k)];
  iota = diag([fi1,fi2]);
  uw = -mat.rho{k}*w^2*fi0*[rim(k)^3; rim(k+1)^3];
  G = [rim(k)^kappa rim(k)^-kappa; rim(k+1)^kappa rim(k+1)^-kappa];
  C = (G*iota)\(u - uw);

  % Discretize radius vector
  rv = linspace(rim(k), rim(k+1), rdiv);
  rvstart = (k-1)*rdiv + 1;
  rvend = k*rdiv;
  rArr(rvstart:rvend) = rv;

  % Calculate discrete displacement vector
  dv = -mat.rho{k}*w^2*fi0*rv.^3 + C(1)*fi1*rv.^kappa + C(2)*fi2*rv.^-kappa;
  uArr(b,rvstart:rvend) = dv; % Discrete displacement throughout the rim

  % Strain
  e1 = dv ./ rv;
  e3 = -3*mat.rho{k}*w^2*fi0*rv.^2 + kappa*(C(1)*fi1*rv.^(kappa-1) - C(2)*fi2*rv.^(-kappa-1));
  e2 = zeros(size(e1)); % no strain in the axial or shear directions
  e4 = zeros(size(e1));
  eArr = [e1; e2; e3; e4]; % strain in each direction [hoop, axial, raidal, shear]

  % Stress
  stress = mat.Q{b,k} * eArr;
  sArr(:,rvstart:rvend,b) = stress;
end

if getappdata(prog,'Canceling')
  delete(prog)
  vari = -1;
  return
end
perc = (b / vari);
waitbar(perc,prog)


delete(prog)
end
