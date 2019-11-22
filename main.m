clc
clear
close('all','force')

%% -----------------------------------------------------------------------------
%  Define global variables
%  -----------------------------------------------------------------------------
% Global variables and arrays
global U w rim arraySize rArr uArr sArr eArr rdiv vari
% Global structures
global mat plotWhat

%% -----------------------------------------------------------------------------
%  Define initial conditions and rotor size
%  -----------------------------------------------------------------------------
% simulation type:
  % pe = steady state perfectly elastic
  % ve = steady state viscoelastic

st = 'pe';
Ftype = 'MaxR'; % Options: TsaiWu, MaxR

% Rotor
rim = [.10, .110, 0.17]; % rim radii in [m]
rdiv = 30; % number of points per rim to analyze
delta = [0.4, 0]/1000; % [mm]
sigb = [0, 0];
mats = {'Alumin_7075_t6.mat','Glass_Epoxy_Ha1999.mat'};
compFunc = @IM7_8552_Tzeng2001; % compliance function, input 'no' to turn off creep modeling

% Time
simTime = 1;
timeUnit = 's'; % s = sec, h = hours, d = days
numberOfSteps = 3;
startime = 1;

% Velocity
iRPM = 38000; % Initial rpm
dThicc = 0.0015; % Damaged ring thickness [m]
degStiffPerc = 0.01; % Degraded stiffness percent
vdiv = 1; % number of points to analyze between each fixed velocity
failure = false;

% Plotting
plotWhat.rims = rim;
plotWhat.custom1 = 'no';

plotWhat.radDis = 'no';         % Radial displacement v. radius
plotWhat.radStr = 'no';         % Radial stress v. radius
plotWhat.hoopStr = 'no';        % Hoop stress v. radius
plotWhat.strengthRatio = 'no';   % Strength Ratio v. radius

plotWhat.disGif = 'no';          % Displacement gif, surface plot
plotWhat.disGifName = 'Displacement.gif';
plotWhat.radGif = 'no';          % Radial stress gif, surface plot
plotWhat.radialGifName = 'Radial Stress.gif';
plotWhat.hoopGif = 'no';         % Hoop stress gif, surface plot
plotWhat.hoopGifName = 'Hoop Stress.gif';
plotWhat.interval = 1;           % Display time interval on figures
plotWhat.delay = 0;              % Time delay in seconds between frames in the gifs,
                                 %   0 is fastest

%% -----------------------------------------------------------------------------
%  Start Program
%  -----------------------------------------------------------------------------
fprintf('Material behavior: %s\n',st)
fprintf('Failure model: %s\n', Ftype)
fprintf('Simulation time: %1.0f %s\n',simTime,timeUnit)
fprintf('Number of rims: %1.0f\n',length(rim)-1)
fprintf('Material selections: %s\n', mats{1})
fprintf('                     %s\n', mats{2:end})
fprintf('Initial rotational velocity: %1.0f rpm\n\n', iRPM)
fprintf('Program Start: \n')

%% -----------------------------------------------------------------------------
%  Check input variables
%  -----------------------------------------------------------------------------
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

  if length(iRPM) > 1
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
fprintf('Create Variable Arrays: Complete\n')
%% -----------------------------------------------------------------------------
%  Preallocate variables
%  -----------------------------------------------------------------------------
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
%  Create Q matrices for all materials
%  -----------------------------------------------------------------------------
b = 1;

for k = 1:length(mats)
  mat.file{k} = ['MaterialProperties\', mats{k}];
  matProp = load(mat.file{k});
  mat.Q{b,k} = stiffMat(matProp.mstiff, compFunc);
  mat.rho{k} = matProp.rho;

  try
    mat.stren{k} = matProp.stren;
  catch
    if ~strcmp(Ftype, 'none')
      error('Yield Strength for one or more materials is not specified. Can not complete the selected simulation.')
    end
    break
  end
end

fprintf('Create Material Property Matrices: Complete\n')

iter = 0; % 0 corresponds to the initial starting time and velocity. This might change to vari when incorporating VE behavior

while ~failure
%% -----------------------------------------------------------------------------
%  Current time and velocity
%  -----------------------------------------------------------------------------
  cRPM = iRPM + 10*iter;
  w = (pi/30) * cRPM;

  fprintf('Faiure = %f\n', failure);
  fprintf('Iteration = %f\n', iter);
  fprintf('cRPM = %f\n\n', cRPM)
%   cont = input('Would you like to contine? ');
%   if ~cont
%     fprintf('Program Ended\n');
%     return
%   end

%% -----------------------------------------------------------------------------
%  Displacement of rim surfaces
%  -----------------------------------------------------------------------------
% Calculate displacement magnitude at the inner and outer surface of each rim
% these are used as boundary conditions to find C. ~ is used to disregard
% output of force vector results. These can be important for debugging and
% verification purposes, but are not necessary for the program. Check function
% discription for mor info
  [~, ~, ~, ~] = boundaryConditions(sigb, delta);

  fprintf('Calculate Boundary Conditions: Complete\n')
%% -----------------------------------------------------------------------------
%  Rotor stress strain calculations
%  -----------------------------------------------------------------------------
% Calculate discrete displacement, stain, and stress for each rim ~ here is
% used to the [C] matrix output. This is useful for debugging and
% verification purposes but not necessary for the function. Check function
% description for mor info
  [~] = discretizeStressStrain(delta);

  fprintf('Descretize Stress/Strain: Complete\n')
%% -----------------------------------------------------------------------------
%  Failure behavior and locations
%  -----------------------------------------------------------------------------
% Failure index and type calculations. Calcuates the failure index using the
% selected faliure modes, determines the type of failure (cracking or burst),
% and identifies the failure location. Outputs determine if the simulation
% should continue or end.

  if ~strcmp(Ftype, 'none')
    [failure, ~, Fmode, Floc] = failureIndex(Ftype);
    fprintf('Failure Analysis: Complete\n');
    [newRotor, rimInd] = degradeRotor(Floc, dThicc);

    mat.file = ['MaterialProperties\', mats{rimInd}];
    ringProps = load(mat.file);
    ringProps.mstiff(2) = degStiffPerc * matProp.mstiff(2);
    ringStiff = stiffMat(ringProps.mstiff, compFunc);
    ringRho = ringProps.rho;
    ringStren = rinProps.stren;

    tempQ = cell(1,length(newRotor));
    tempQ{1:rimInd+2} = {mat.Q{1:rimInd}, ringStiff, mat.Q{rimInd}};

    if
    tempQ{rimInd+2} = {mat.Q{rimInd}};

    %Integrate ring stiffness, density, and strength into the larger mat structure.
  end

  iter = iter + 1;
end
%% -----------------------------------------------------------------------------
%  Make Plots
%  -----------------------------------------------------------------------------
  plotStressStrain()

fprintf('Create Output Plots: Complete\n\n')
fprintf('Program Complete\n')
