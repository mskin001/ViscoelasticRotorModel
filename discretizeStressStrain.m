function [C] = discretizeStressStrain()
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
global rotor numRims delta U w rArr uArr sArr eArr rdiv

%% -----------------------------------------------------------------------------
% Calculate stress, strain, displacement for each rim
% ------------------------------------------------------------------------------
for k = 1:numRims
  rim = rotor.radii{k};
  Q11 = rotor.Q{k}(1,1);
  Q13 = rotor.Q{k}(1,3);
  Q33 = rotor.Q{k}(3,3);
  kappa = sqrt(Q11/Q33); % intermediate variable of stiffness ratio

  % intermediate variables for calculating the stiffness matrix
  fi0 = 1/(Q33 * (9 - kappa^2));
  fi1 = 1/(Q13 + (kappa * Q33));
  fi2 = 1/(Q13 - (kappa * Q33));

  % Calculate absolute displacement at the inner and outer surface of each rim
  u = [U(1,k); U(1,k+1) - delta(k)];
  iota = diag([fi1,fi2]);
  uw = -rotor.rho(k)*w^2*fi0*[rim(1)^3; rim(2)^3];
  G = [rim(1)^kappa rim(1)^-kappa; rim(2)^kappa rim(2)^-kappa];
  C = (G*iota)\(u - uw);

  % Discretize radius vector
  rv = linspace(rim(1), rim(2), rdiv);
  rvstart = (k-1)*rdiv + 1;
  rvend = k*rdiv;
  rArr(rvstart:rvend) = rv;

  % Calculate discrete displacement vector
  dv = -rotor.rho(k)*w^2*fi0*rv.^3 + C(1)*fi1*rv.^kappa + C(2)*fi2*rv.^-kappa;
  uArr(rvstart:rvend) = dv; % Discrete displacement throughout the rim

  % Strain
  e1 = dv ./ rv;
  e3 = -3*rotor.rho(k)*w^2*fi0*rv.^2 + kappa*(C(1)*fi1*rv.^(kappa-1) - C(2)*fi2*rv.^(-kappa-1));
  e2 = zeros(size(e1)); % no strain in the axial or shear directions
  e4 = zeros(size(e1));
  eArr = [e1; e2; e3; e4]; % strain in each direction [hoop, axial, raidal, shear]

  % Stress
  stress = rotor.Q{k} * eArr;
  sArr(:,rvstart:rvend) = stress;
end

end
