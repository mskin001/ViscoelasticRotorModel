function [outerRadii] = degradeRotor(dThicc)
% This function takes in the failure location and adds in a ring of degraded
% material with a thickness of dThicc and degraded material properties
% degMatProp. This is accomplished by identifying the rim where the failure
% occurs, dividing that into two smaller rims separated by a degraded ring of
% material with inner radius degInner and outer radius degOuter, and adding
% those radii to the overall rotor dimensions (rim).

global rotor failure delta

innerPos = failure.rim; % position of the failed rim in the rotor
ringPos = innerPos + 1;
outerPos = innerPos + 2;

ringInner = failure.loc - dThicc/2; % Inner radius of degraded ring
ringOuter = failure.loc + dThicc/2; % Outer radius of degraded ring

% Dimension and position of each new rim made from segmenting the origional
innerRadii = [rotor.radii{innerPos}(1), ringInner];
ringRadii = [ringInner, ringOuter]; % Degrated ring radii. This vector will always be size [1,2]
outerRadii = [ringOuter, rotor.radii{innerPos}(2)];

rotor.radii(outerPos:end+2) = rotor.radii(innerPos:end);
rotor.radii{innerPos} = innerRadii;
rotor.radii{ringPos} = ringRadii;
rotor.radii{outerPos} = outerRadii;

delta(outerPos:end+2) = delta(innerPos:end);
delta(innerPos:outerPos) = 0;

% Lable intact and failed rims
rotor.intact(outerPos:end+2) = rotor.intact(innerPos:end);
rotor.intact(innerPos) = 1;
rotor.intact(ringPos) = 0;

rotor.matInd(outerPos:end+2) = rotor.matInd(innerPos:end);
rotor.matInd(outerPos) = rotor.matInd(innerPos);
rotor.matInd(ringPos) = rotor.matInd(outerPos);
