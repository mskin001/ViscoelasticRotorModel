function [Fw, Fe0, Fe1, Fd, Fb, K] = boundaryConditions(E, sigb, delta)
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
global U rim w mat arraySize vari

K = zeros(arraySize,arraySize);    % global stiffness matrix
Fw = zeros(arraySize,1);           % vector of centrifugal forces
Fe0 = zeros(arraySize,1);          % vector of axial strain forces
Fe1 = zeros(arraySize,1);          % vector of axial strain forces
Fd = zeros(arraySize,1);           % vector of presress forces
Fb = zeros(arraySize,1);           % vector of boundary condition forces
%% -----------------------------------------------------------------------------
% Calculate boundary displacement at the inner and outer surface of each rim
% ------------------------------------------------------------------------------
prog = waitbar(0,'Calculating Boundary Conditions', 'CreateCancelBtn',...
  'setappdata(gcbf,''Canceling'',1)');
setappdata(prog,'Canceling',0);

if length(w) > 1 % First part is for qdve simulation
  % Under construction
elseif length(w) == 1 % Second part for pe and ve simulations
  for b = 1:vari
    for k = 1:length(rim)-1
      [Q, kappa, fi] = findMatPropConsts(b,k);

      z = rim(k)/rim(k+1);
      z1 = z^-kappa + z^kappa;
      z2 = z^-kappa - z^kappa;
      
      % Local stiffness matrix
      kMat = (1/z2) * [kappa*z1*Q(3,3)-z2*Q(1,3), -2*kappa*Q(3,3);
                      -2*kappa*Q(3,3), kappa*z1*Q(3,3)+z2*Q(1,3)];

      % Global stiffness matrix for the entire
      K(k:k+1, k:k+1) = K(k:k+1, k:k+1) + kMat;
      
      % Calculate force terms
      uw = -(mat.rho{k})*w^2*fi(1)*[rim(k)^3; rim(k+1)^3];
      u0 = fi(5)*[rim(k); rim(k+1)];
      u1 = fi(4)*[rim(k)^2; rim(k+1)^2];
      
      fsig = -(mat.rho{k})*w^2*fi(6)*[-rim(k)^3; rim(k+1)^3];
      f0 = fi(8)*E(b,1)*[-rim(k); rim(k+1)];
      f1 = fi(7)*E(b,2)*[-rim(k)^2; rim(k+1)^2];
            
      fw = -fsig + kMat*uw;
      fe0 = -f0 + kMat*u0;
      fe1 = -f1 + kMat*u1;
      
      Fw(k:k+1) = Fw(k:k+1) + fw;
      Fe0(k:k+1) = Fe0(k:k+1) + fe0;
      Fe1(k:k+1) = Fe1(k:k+1) + fe1;
      Fd(k:k+1) = Fd(k:k+1) + kMat*[0;delta(k)];
    end

    % Force from external pressure applied to inner and outer surface of the rim
    Fb(1) = -rim(1)*sigb(1);
    Fb(end) = rim(end)*sigb(end);
    % Displacement at the inner and outer radius of each rim
    U(b,:) = K \ (Fb + Fw + Fe0 + Fe1 + Fd);

    % Reset matrix values
    K = K .* 0;
    Fw = Fw .* 0;
    Fb = Fb .* 0;
    Fd = Fd .* 0;
    
    if getappdata(prog,'Canceling')
      delete(prog)
      vari = -1;
      return
    end
    perc = (b / vari);
    waitbar(perc,prog) 
  end
  
else
  error('Could not identify intented simulation. Review input variables in "main.m"')
end

delete(prog)
end
