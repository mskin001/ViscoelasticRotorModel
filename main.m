% clc
clear
close('all','force')
%% -----------------------------------------------------------------------------
% Define global variables
% ------------------------------------------------------------------------------
% Global variables and arrays
global U w rim arraySize rArr uArr sArr eArr t vari
% Global structures
global mat plotWhat

%% -----------------------------------------------------------------------------
% Define initial conditions and rotor size
% ------------------------------------------------------------------------------
% simulation type:
  % pe = steady state perfectly elastic
  % ve = steady state viscoelastic
  % qdve = quasi-dynamic viscoelasticu
st = 'pe';

% Rotor
% rim = [0.05; 0.1]; % single rim Ha 1999
rim = [0.03789, 0.079];
% rim = [0.03786, 0.08393, 0.14707];
% rim = [0.0762, .1524]; % Tzeng2001
rdiv = 30; % number of points per rim to analyze
delta = [0]/1000; % [mm]
sigb = [-.778e8, 0];
mats = {'Glass_Epoxy_Ha1999.mat'};
% mats = {'T300_2500_Ha2001.mat'};
% mats = {'AS_H3501_Ha1999.mat'; 'IM6_Epoxy_Ha1999.mat'};

% Time/creep
simTime = 10;
timeUnit = 's'; % s = sec, h = hours, d = days
numberOfSteps = 3;
compFunc = @IM7_8552_Tzeng2001; % compliance function, input 'no' to turn off creep modeling

% Speed/velocity
rpm = 60000;
vdiv = 1; % number of points to analyze between each fixed velocity

% Plotting
plotWhat.rims = rim;
plotWhat.custom1 = 'yes';

plotWhat.disGif = 'no';          % Displacement gif, surface plot
plotWhat.disGifName = 'Displacement.gif';
plotWhat.radDis = 'no';

plotWhat.radGif = 'no';          % Radial stress gif, surface plot
plotWhat.radialGifName = 'Radial Stress.gif';
plotWhat.radStr = 'yes';         % Radial stress v. radius plot

plotWhat.hoopGif = 'no';         % Hoop stress gif, surface plot
plotWhat.hoopGifName = 'Hoop Stress.gif';
plotWhat.hoopStr = 'yes';        % Hoop stress v. radius plot

plotWhat.axialStr = 'no';

plotWhat.interval = 1;           % Display time interval on figures
plotWhat.delay = 0;              % Time delay in seconds between frames in the gifs,
                                 %   0 is fastest
%% -----------------------------------------------------------------------------
% Start Program
% ------------------------------------------------------------------------------
fprintf('Simulation type: %s\n',st)
fprintf('Simulation time: %1.0f %s\n',simTime,timeUnit)
fprintf('Number of rims: %1.0f\n',length(rim)-1)
fprintf('Material Selections: %s\n', mats{1:end})
% fprintf('                     %s\n', mats{2:end})
fprintf('Rotational velocity: %1.0f rpm\n\n', rpm)
fprintf('Program Start: \n')

%% -----------------------------------------------------------------------------
% Check input variables
% ------------------------------------------------------------------------------
% General
if length(rim)-1 ~= length(mats) || length(rim)-1 ~= length(delta)
  error('Error in rim, mat, and delta. There must be an equal number of rims, materials, and interferance values.\n')
end

for k = 1:length(mats)
    m = ['MaterialProperties\', mats{k}];
    mp = load(m);
    try
        verified = mp.verified;
        if ~strcmp(verified, 'True')
            warning('The material %s has not been verified\n', mats{k})
        end
    catch
        warning('The material %s has not been verified\n', mats{k})
    end
end

% Simulation specific
if strcmp(st,'pe')
  simTime = 1; % steady state = not time changes
  compFunc = 'no'; % Redefineds compFunc to reflect a constant elastic matrix

  if length(rpm) ~= 1
    error('Rotational velocity not specified. Please specify a single veloctiy\n')
  end

elseif strcmp(st,'ve')
  if simTime < 2
    error('Simulation time is less than 2 units. Please specify a simulation time larger than 2, or change to a perfectly elastic simulation type by defining st = ''pe'' \n')
  end

  if length(rpm) > 1
    error('Rotational velocity not specified. Please specify a single veloctiy\n')
  end

  if ~contains(['s','h','d'], timeUnit)
    warning('This unit of time is not supported, supported units are to seconds, hours, or days\n')
    disp('Ignore the above warning if intentionally using a different unit, check compliance funtion inputs and plotting outputs\n')
  end

  if length(compFunc) ~= length(mats)
    warning('More materials are specified than compliance functions.\n')
    disp('Ignore the above warning if 2 or more materials are the same, or if 1 or more material elasticity is constant\n')
  end

elseif strcmp(st,'qdve')
  error('Simulation type not supported\n')

  if length(rpm) > 1 && vdiv < 2 %#ok<UNRCH>
    warning('vdiv is less than 2. To study changes in velocity this should be 2 or larger\n')
  end
else
  error('Simulation type not specified or not supported.\n')
end

fprintf('Check Input Variables: Complete\n')
%% -----------------------------------------------------------------------------
% Create speed/time arrays depending on simulation global
% ------------------------------------------------------------------------------
if length(rpm) > 1
  wComp = zeros(vdiv, 1);
  for k = 1:length(rpm)-1
    w1 = (pi/30) * rpm(k);
    w2 = (pi/30) * rpm(k+1);
    wStart = (k-1)*vdiv + 1;
    wEnd = k*vdiv;
    wComp(wStart:wEnd) = linspace(w1, w2, vdiv);
  end
  w = wComp;
  vari = length(wComp);
elseif simTime > 1
  if strcmp(timeUnit, 'h')
    simTime = simTime * 3600; % convert hours to seconds
  elseif strcmp(timeUnit, 'd')
    simTime = simTime * 24 * 3600; % Convert days to seconds
  end
  tArr = [1, 1*10e5, 1*10e10]; % Assumes 1 sec time intervals
  w = (pi/30) * rpm;
  vari = length(tArr);
  addpath('ComplianceFunctions')
else
  w = (pi/30) * rpm;
  vari = 1;
  tArr = 1;
end

fprintf('Create Variable Arrays: Complete\n')
%% -----------------------------------------------------------------------------
% Preallocate variables
% ------------------------------------------------------------------------------
arraySize = length(rim);
U = zeros(vari,arraySize);
rArr = zeros(1,(arraySize-1)*rdiv);    % radius vector for descretization
uArr = zeros(vari,(arraySize-1)*rdiv);    % displacement vector for discretization
sArr = zeros(4,(arraySize-1)*rdiv,vari);    % stress vector
eArr = zeros(4, rdiv);    % strain vector in each direction
mat = struct();
mat.file = cell(length(mats),1);
mat.Q = cell(vari,length(mats));

fprintf('Preallocate Memory: Complete\n')
%% -----------------------------------------------------------------------------
% Create Q matrices for all materials
% ------------------------------------------------------------------------------
prog = waitbar(0,'Creating Material Property Matrices', 'CreateCancelBtn',...
  'setappdata(gcbf,''Canceling'',1)');
setappdata(prog,'Canceling',0);

for b = 1:vari
  t = tArr(b);  
  if getappdata(prog,'Canceling')
    delete(prog)
    return
  end
  
  for k = 1:length(mats)
      mat.file{k} = ['MaterialProperties\', mats{k}];
      matProp = load(mat.file{k});
      mat.Q{b,k} = stiffMat(matProp.mstiff, compFunc);
      mat.rho{k} = matProp.rho;
      
      try
        mat.stren{k} = matProp.stren;
      catch
        break
      end
  end
  
  perc = (b / vari);
  waitbar(perc,prog)  
end

delete(prog)
fprintf('Create Material Property Matrices: Complete\n')

%% Find Constants
[E,~] = findAxialStrainCoeff(sigb);

%% ----------------------------------------------------------------------------
% Calculate displacement magnitude at the inner and outer surface of each rim
% these are used as boundary conditions to find C. ~ is used to disregard
% output of force vector results. These can be important for debugging and
% verification purposes, but are not necessary for the program. Check function
% discription for mor info
[~, ~, ~, ~, ~, ~] = boundaryConditions(E,sigb, delta);

if vari == -1
  return
end
fprintf('Calculate Boundary Conditions: Complete\n')

%% -----------------------------------------------------------------------------
% Calculate discrete displacement, stain, and stress for each rim ~ here is
% used to the [C] matrix output. This is useful for debugging and
% verification purposes but not necessary for the function. Check function
% description for mor info
[~] = discretizeStressStrain(rdiv, delta, E, sigb);

if vari == -1
  return
end
fprintf('Descretize Stress/Strain: Complete\n')
%% -----------------------------------------------------------------------------
% Make Plots
plotStressStrain()

fprintf('Create Output Plots: Complete\n\n')
fprintf('Program Complete\n')
