function [C] = discretizeStressStrain(rdiv, delta, C1, C2, E0, E1)
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
  % Under construction
elseif length(w) == 1
  for b = 1:vari
    for k = 1:length(rim)-1
      [~, kappa, fi] = findMatPropConsts(1,1);

      % Discretize radius vector
      rv = linspace(rim(k), rim(k+1), rdiv);
      rvstart = (k-1)*rdiv + 1;
      rvend = k*rdiv;
      rArr(rvstart:rvend) = rv;

      % Calculate discrete displacement vector
      C = [C1,C2];
      dv = -mat.rho{k}*w^2*fi(1)*rv.^3 + C(1)*fi(2)*rv.^kappa + C(2)*fi(3)*rv.^-kappa ...
            + fi(4)*E1*rv.^2 + fi(5)*E0*rv;
      uArr(b,rvstart:rvend) = dv; % Discrete displacement throughout the rim

      % Strain
      e1 = dv ./ rv;
      e3 = -3*mat.rho{k}*w^2*fi(1)*rv.^2 + kappa*C(1)*fi(2)*rv.^(kappa-1) - kappa*C(2)*fi(3)*rv.^(-kappa-1) ...
              + 2*fi(4)*E1*rv + fi(5)*E0;
      e2 = E0 + E1*rv; % no strain in the axial or shear directions
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
  end
end

delete(prog)
end
