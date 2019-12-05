clc
clear
close('all','force')
warning('off','all')
%% -----------------------------------------------------------------------------
%  Declare global variables
%  -----------------------------------------------------------------------------
% Global variables and arrays
global U w numRims arraySize rArr uArr sArr eArr rdiv
% Global structures
global plotWhat rotor

%% -----------------------------------------------------------------------------
%  Define initial conditions and rotor size
%  -----------------------------------------------------------------------------
% simulation type:
  % pe = steady state perfectly elastic
  % ve = steady state viscoelastic
st = 'pe';
Ftype = 'TsaiWu'; % Options: TsaiWu, MaxR

% Rotor
rim = [0.001, .055, 0.080]; % rim radii in [m]
rdiv = 60; % Number of points per rim to analyze
delta = [0.175, 0]/1000; % [mm]
sigb = [-67.5e6, 0];
mats = {'Al7075T651_Corbin2005.mat','G30-500_8604_Corbin2005.mat'};
compFunc = @IM7_8552_Tzeng2001; % Compliance function, input 'no' to turn off creep modeling
dThicc = 0.0015; % Damaged ring thickness [m]
degStiffPerc = 0.005; % Stiffness degraded percent
failure = false; % No initial failures
Fmode = 'none'; % No initial failure mode
burstSpeed = [];
matrixFailureLoc = [];

% Time
simTime = 1;
timeUnit = 's'; % s = sec, h = hours, d = days
numberOfSteps = 3; % Number of time steps
startime = 1;

% Velocity
iRPM = 65000; % Initial rpm
vdiv = 1; % number of points to analyze between each fixed velocity


% Plotting
plotWhat.rims = rim;
plotWhat.custom1 = 'no';

plotWhat.radDis = 'yes';         % Radial displacement v. radius
plotWhat.radStr = 'yes';         % Radial stress v. radius
plotWhat.hoopStr = 'yes';        % Hoop stress v. radius
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

%% -----------------------------------------------------------------------------
%  Preallocate variables
%  -----------------------------------------------------------------------------
numRims = length(rim)-1;
rotor = struct();
rotor.matInd = linspace(1, numRims, numRims);
rotor.radii = cell(1, numRims);
rotor.Q = cell(1, numRims);
rotor.rho = zeros(1,numRims);
rotor.stren = cell(1,numRims);
rotor.intact = ones(1,numRims); % '1' indicates the rim is intact. All rims initially intact.

fprintf('Preallocate Memory: Complete\n')

for k = 1:numRims
  rotor.radii{k} = [rim(k),rim(k+1)];
end

iter = 0; % 0 corresponds to the initial starting time and velocity. This might change to vari when incorporating VE behavior
while ~strcmp('Burst',Fmode)
  arraySize = numRims+1; % Intermediate variable
  U = zeros(1,arraySize);
  rArr = zeros(1,(arraySize-1)*rdiv);    % radius vector for descretization
  uArr = zeros(1,(arraySize-1)*rdiv);    % displacement vector for discretization
  sArr = zeros(4,(arraySize-1)*rdiv);    % stress vector
  eArr = zeros(4, rdiv);    % strain vector in each direction

  while ~failure
  % -----------------------------------------------------------------------------
  %  Current time and velocity
  %  ---------------------------------------------------------------------------
    cRPM = iRPM + 10*iter;
    w = (pi/30) * cRPM;

    %% -------------------------------------------------------------------------
    %  Create Q matrices and proprogate rotor structures
    %  -------------------------------------------------------------------------
    for k = 1:numRims
      matFile = ['MaterialProperties\', mats{rotor.matInd(k)}];
      matProp = load(matFile);
      rotor.rho(k) = matProp.rho;

      if rotor.intact(k)
        rotor.Q{k} = stiffMat(matProp.mstiff, compFunc);
      else
        % Do something for broken rims
        matProp.mstiff(2) = degStiffPerc * matProp.mstiff(2);
        rotor.Q{k} = stiffMat(matProp.mstiff, 'no');
      end

      try
        rotor.stren{k} = matProp.stren;
      catch
        if ~strcmp(Ftype, 'none')
          error('Yield Strength for one or more materials is not specified. Can not complete the selected simulation.')
        else
        warning('Yield strength not specified for this material')
        fprintf('%s\n', mats{k});
        end
      end

    end

  %% -----------------------------------------------------------------------------
  %  Displacement of rim surfaces
  %  -----------------------------------------------------------------------------
  % Calculate displacement magnitude at the inner and outer surface of each rim
  % these are used as boundary conditions to find C. ~ is used to disregard
  % output of force vector results. These can be important for debugging and
  % verification purposes, but are not necessary for the program. Check function
  % discription for mor info
    [~, ~, ~, ~] = boundaryConditions(sigb, delta);
%     fprintf('Calculate Boundary Conditions: Complete\n')
  %% -----------------------------------------------------------------------------
  %  Rotor stress strain calculations
  %  -----------------------------------------------------------------------------
  % Calculate discrete displacement, stain, and stress for each rim ~ here is
  % used to the [C] matrix output. This is useful for debugging and
  % verification purposes but not necessary for the function. Check function
  % description for mor info
    [~] = discretizeStressStrain(delta);
%     fprintf('Descretize Stress/Strain: Complete\n')
  %% -----------------------------------------------------------------------------
  %  Failure behavior and locations
  %  -----------------------------------------------------------------------------
  % Failure index and type calculations. Calcuates the failure index using the
  % selected faliure modes, determines the type of failure (cracking or burst),
  % and identifies the failure location. Outputs determine if the simulation
  % should continue or end.    
    if ~strcmp(Ftype, 'none')
      [failure, strengthRatio, Fmode, Floc] = failureIndex(Ftype);
%       fprintf('Failure Analysis: Complete\n');
    else
      failure = 1; % When not evaluating failure, set to 1 to iterate exactly once.
      Fmode = 'Burst';
      Floc = 0;
    end

    iter = iter + 1;
    results.rotor{iter} = rotor;
    results.rArr{iter} = rArr;
    results.uArr{iter} = uArr;
    results.sArr{iter} = sArr;
  end

  if ~strcmp('Burst',Fmode)
    [delta] = degradeRotor(Floc, dThicc, delta);
  end
  
  burstSpeed(end+1) = cRPM;
  matrixFailureLoc(end+1) = Floc;
  failure = 0;
  numRims = length(rotor.radii);
end


%% -----------------------------------------------------------------------------
%  Make Plots
%  -----------------------------------------------------------------------------
plotStressStrain()

fprintf('Create Output Plots: Complete\n\n')
fprintf('Program Complete\n')
