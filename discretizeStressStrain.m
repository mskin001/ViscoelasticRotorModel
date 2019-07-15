function [C] = discretizeStressStrain(rdiv, delta, e0, e1)
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

if length(w) > 1
% Under Construction
elseif length(w) == 1
  for b = 1:vari
    for k = 1:length(rim)-1
      [~, kappa, fi] = findMatPropConsts(b,k);
      
      % Calculate absolute displacement at the inner and outer surface of each rim
      u = [U(b,k); U(b,k+1) - delta(k)];
      iota = diag([fi(2),fi(3)]);
      uw = -mat.rho{k}*w^2*fi(1)*[rim(k)^3; rim(k+1)^3];
      u0 = fi(5)*[rim(k); rim(k+1)];
      u1 = fi(4)*[rim(k)^2; rim(k+1)];
      G = [rim(k)^kappa rim(k)^-kappa; rim(k+1)^kappa rim(k+1)^-kappa];
      C = (G*iota)\(u - uw - u0 - u1);

      % Discretize radius vector
      rv = linspace(rim(k), rim(k+1), rdiv);
      rvstart = (k-1)*rdiv + 1;
      rvend = k*rdiv;
      rArr(rvstart:rvend) = rv;

      % Calculate discrete displacement vector
      dv = -mat.rho{k}*w^2*fi(1)*rv.^3 + C(1)*fi(2)*rv.^kappa + C(2)*fi(3)*rv.^-kappa...
           +fi(4)*e1*rv.^2 + fi(5)*e0*rv;
      uArr(b,rvstart:rvend) = dv; % Discrete displacement throughout the rim

      % Strain
      ep1 = dv ./ rv;
      ep2 = e0 + e1*rv;
      ep3 = -3*mat.rho{k}*w^2*fi(1)*rv.^2 + kappa*(C(1)*fi(2)*rv.^(kappa-1) - C(2)*fi(3)*rv.^(-kappa-1));
      ep4 = zeros(size(ep1));
      eArr = [ep1; ep2; ep3; ep4]; % strain in each direction [hoop, axial, raidal, shear]

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
  end
end

delete(prog)
end
