function [Fw, Fd, Fb, K] = boundaryConditions(sigb, delta)
% This funciton calculates the displacement at the inner and outer surface of
% each rim. These will be used as boundary conditions to calculate the stress
% displacement throughout the entire rotor. The function is separated into
% two similar parts to reduce memory requirements. The first part simulates
% steady state simulations at various velocities. The second part simulates
% constant velocity at non constant time, i.e. creep behavior.
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
global rotor numRims U w arraySize

K = zeros(arraySize,arraySize);    % global stiffness matrix
Fw = zeros(arraySize,1);           % vector of centrifugal forces
Fd = zeros(arraySize,1);           % vector of presress forces
Fb = zeros(arraySize,1);           % vector of boundary condition forces

%% -----------------------------------------------------------------------------
% Calculate boundary displacement at the inner and outer surface of each rim
% ------------------------------------------------------------------------------
for k = 1:numRims
  rim = rotor.radii{k};
  Q11 = rotor.Q{k}(1,1);
  Q13 = rotor.Q{k}(1,3);
  Q33 = rotor.Q{k}(3,3);
  kappa = sqrt(Q11/Q33);

  % intermediate variables for calculating the stiffness matrix
  fi0 = 1/(Q33 * (9 - kappa^2));
  fi3 = (3 * Q33 + Q13) / (Q33 * (9 - kappa^2));

  z = rim(1)/rim(2);
  z1 = z^-kappa + z^kappa;
  z2 = z^-kappa - z^kappa;

  % Local stiffness matrix
  kMat = (1/z2) * [kappa*z1*Q33-z2*Q13, -2*kappa*Q33;
                  -2*kappa*Q33, kappa*z1*Q33+z2*Q13];

  % Global stiffness matrix for the entire
  K(k:k+1, k:k+1) = K(k:k+1, k:k+1) + kMat;

  fsig = -(rotor.rho(k))*w^2*fi3*[-rim(1)^3; rim(2)^3];
  uw = -(rotor.rho(k))*w^2*fi0*[rim(1)^3; rim(2)^3];

  fw = -fsig + kMat*uw;
  Fw(k:k+1) = Fw(k:k+1) + fw;

  Fd(k:k+1) = Fd(k:k+1) + kMat*[0;delta(k)];
end

% Force from external pressure applied to inner and outer surface of the rim
Fb(1) = -rotor.radii{1}(1)*sigb(1);
Fb(end) = rotor.radii{end}(2)*sigb(end);
% Displacement at the inner and outer radius of each rim
U(:) = K \ (Fb + Fw + Fd);

% Reset matrix values
K = K .* 0;
Fw = Fw .* 0;
Fb = Fb .* 0;
Fd = Fd .* 0;

end
