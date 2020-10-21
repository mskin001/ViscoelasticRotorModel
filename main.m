clc
clear all
close('all','force')
format long
%% -----------------------------------------------------------------------------
% Define global variables
% ------------------------------------------------------------------------------
% Global variables and arrays
global U w rim arraySize rArr uArr sArr eArr tauArr vari
% Global structures
global mat plotWhat results

%% -----------------------------------------------------------------------------
% Define initial conditions and rotor size
% ------------------------------------------------------------------------------
% simulation type:
  % pe = steady state perfectly elastic
  % ve = steady state viscoelastic
st = 'pe';

% Rotor
% rim = [0.03789; 0.07901]; % single rim Ha 1999
% rim = [.1, 0.8];
rim = [.05, .1];
% rim = [0.08, 0.2]; % Perez-Aparicio 2011
% rim = [0.0762, .1524]; % Tzeng2001
rdiv = 30; % number of points per rim to analyze
delta = [0]/1000; % [mm]
sigb = [0, 0];
% mats = {'salehian_Incl718.mat'};
mats = {'CFRP_Aparicio2011.mat'};

% Time/creep
tmax = 20; %seconds?
tStep = 1; %second between steps
simTime = tmax;
timeUnit = 's'; % s = sec, h = hours, d = days
compFunc = {'no' 'no'}; % compliance function, input 'no' to turn off creep modeling

% Speed/velocity
rpm = 50000;
% rpmMax = 97682;
accType = 'const';

% The following if statment controls which acceleration function is used
% based on the variable string accType. This variable is also referenced in
% main.m and shearStress.m. Apply changes with caution.
if strcmp(accType, 'const')
  % Constant
  alpha = @(t,wIni) 250 * t; % gets multiplied by tStep at end of while loop
  % alpha = @(t,wIni) 3.6e6 * t;
elseif strcmp(accType, 'Linear')
  % Linear Acceleration:
  alpha = @(t,wIni) wIni + 230*t; % gets multiplied by b*tStep at end of while loop
elseif strcmp(accType, 'Exponential')
  % Exponential growth: (Use this one)
  alpha = @(t,wIni) wIni + 1.089^(5*t);
else
  % Sin behavior:
  alpha = @(t,wIni) sin: wIni + 2356.2*sin((2*pi/40)*t + (3*pi/2)) + 2356.2;
end
%-------------------
% %Exponential behavior:
% T = tmax / log(rpmMax/rpm);
% alpha = @(t,wIni) (wIni * exp(t/T));

% Plotting
legTxt = {'Current model', 'Aparicio 2011'};
% legTxt = {'0 sec', '5 sec', '10 sec', '15 sec', '20 sec'}; % Controls legend entries for graphs
plotWhat.custom1 = 'no';        % any custom plot. Go to plotStressStrain.m to modify (first if statement)
plotWhat.radDis = 'no';          % Radial displacement v. radius
plotWhat.radStr = 'no';         % Radial stress v. radius plot
plotWhat.hoopStr = 'no';        % Hoop stress v. radius plot
plotWhat.shearStr = 'yes';       % Shear stress v. radius
plotWhat.peakStr = 'no';        % 2-yaxis plot. Peak stress location and SR v. time

plotWhat.disGif = 'no';          % Displacement gif, surface plot
plotWhat.disGifName = 'Displacement.gif';
plotWhat.radGif = 'no';          % Radial stress gif, surface plot
plotWhat.radialGifName = 'Radial Stress.gif';
plotWhat.hoopGifName = 'Hoop Stress.gif';
plotWhat.hoopGif = 'no';         % Hoop stress gif, surface plot

plotWhat.interval = 25;          % Display time interval on figures
plotWhat.delay = 0;              % Time delay in seconds between frames in the gifs,
                                 %   0 is fastest


%% -----------------------------------------------------------------------------
% Start Program
% ------------------------------------------------------------------------------
fprintf('Simulation type: %s\n',st)
fprintf('Simulation time: %1.0f %s\n',simTime,timeUnit)
fprintf('Number of rims: %2.0f\n',length(rim)-1)
fprintf('Material Selections: %s\n', mats{1:end})
% fprintf('                     %s\n', mats{2:end})
fprintf('Rotational velocity: %6.3f rpm\n\n', rpm)
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
        if ~verified
            warning('The material %s has not been verified\n', mats{k})
        end
    catch
        warning('The material %s has not been verified\n', mats{k})
    end
end

% Simulation specific
if strcmp(st,'pe')
  simTime = 1; % steady state = not time changes
  compFunc = cell(1,length(mats));
  compFunc(1:end) = {'no'}; % Redefineds compFunc to reflect a constant elastic matrix

  if length(rpm) > 1
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

else
  error('Simulation type not specified or not supported.\n')
end

fprintf('Check Input Variables: Complete\n')

%% -----------------------------------------------------------------------------
% Program Begin
% ------------------------------------------------------------------------------
% vel = zeros(1,tmax/tStep);
w = (pi/30) * rpm; %initial angular velocity
w0 = w;
vari = cast(tmax/tStep,'single');
b = 1;

while b*tStep <= tmax
%   fprintf('Create Variable Arrays: Complete\n')
  %% ---------------------------------------------------------------------------
  % Preallocate variables
  % ----------------------------------------------------------------------------
  arraySize = length(rim);
  U = zeros(1,arraySize);
  rArr = zeros(1,(arraySize-1)*rdiv);    % radius vector for descretization
  uArr = zeros(1,(arraySize-1)*rdiv);    % displacement vector for discretization
  sArr = zeros(4,(arraySize-1)*rdiv);    % stress vector
  tauArr = zeros(1,(arraySize-1)*rdiv);
  eArr = zeros(4, rdiv);    % strain vector in each direction

%   fprintf('Preallocate Memory: Complete\n')

  %% ---------------------------------------------------------------------------
  % Create Q matrices for all materials
  % ----------------------------------------------------------------------------

  for k = 1:length(mats)
    func = compFunc{k};
    mat.file{k} = ['MaterialProperties\', mats{k}];
    matProp = load(mat.file{k});
    mat.Q{1,k} = stiffMat(matProp.mstiff, func);
    mat.rho{k} = matProp.rho;
%     mat.rho{k} = @(r) 7800 + 10.*r + 100.*r.^2 + 1000.*r.^3;

    try
      mat.stren{k} = matProp.stren;
    catch
      break
    end
  end

%   fprintf('Create Material Property Matrices: Complete\n')
  %% ---------------------------------------------------------------------------
  % Calculate displacement magnitude at the inner and outer surface of each rim
  % these are used as boundary conditions to find C. ~ is used to disregard
  % output of force vector results. These can be important for debugging and
  % verification purposes, but are not necessary for the program. Check function
  % discription for mor info
  [~, ~, ~, ~] = boundaryConditions(sigb, delta);

%   fprintf('Calculate Boundary Conditions: Complete\n')
  %% ---------------------------------------------------------------------------
  % Calculate discrete displacement, stain, and stress for each rim ~ here is
  % used to the [C] matrix output. This is useful for debugging and
  % verification purposes but not necessary for the function. Check function
  % description for mor info
  [~] = discretizeStressStrain(rdiv, delta);

%   fprintf('Descretize Stress/Strain: Complete\n')

  %%----------------------------------------------------------------------------
  % Calculate the share stress on the rim.
  [~] = shearStress(alpha, accType, b, 0, tStep, rdiv);

  %% ---------------------------------------------------------------------------
  % Store results for post processing
  % ----------------------------------------------------------------------------
  results.uArr{b} = uArr;
  results.sArr{b} = sArr;
  results.tauArr{b} = tauArr;
  results.vel(b) = w; % * (30 / pi);
%   fprintf('Current time: %5.2f\n', b*tStep)
%   fprintf('Iteration %2.0f Complete\n', b)

  %% ---------------------------------------------------------------------------
  % Update angular velocity and time
  % ----------------------------------------------------------------------------

  if strcmp(accType, 'const')
    w = w + alpha(tStep);
  else
    w = alpha(b*tStep,w0);
  end
  results.time(b) = b*tStep;

  b = b + 1;
  
  % ----------------------------------------------------------------------------
  % Calculate rotor energy
  % ----------------------------------------------------------------------------
%   energy(b) = 

end

%% -----------------------------------------------------------------------------
% Calculate failure criterion
% ------------------------------------------------------------------------------
% [SR] = failureIndex(rdiv);

figure()
plot(results.time, results.vel)
hold on
% plot(rArr*1000,SR(15,:), '-o','MarkerIndices', 1:3:length(rArr), 'Linewidth', 1.5)
% plot(rArr*1000,SR(45,:), '--d', 'Color', [0.6350 0.0780 0.1840], 'MarkerIndices', 1:3:length(rArr), 'Linewidth', 1.5)
ylabel('angular velocity')
xlabel('time')
% legend('SR \it t=3', 'SR \it t=9')
% set(gca, 'Fontsize', 12)

%% -----------------------------------------------------------------------------
% Make Plots
% ------------------------------------------------------------------------------
plotStressStrain(legTxt)

% fprintf('Create Output Plots: Complete\n\n')
fprintf('Program Complete\n')
