function [delta] = degradeRotor(rim, Floc, dThicc, degStiffPerc, mat, delta)
% This function takes in the failure location and adds in a ring of degraded
% material with a thickness of dThicc and degraded material properties
% degMatProp. This is accomplished by identifying the rim where the failure
% occurs, dividing that into two smaller rims separated by a degraded ring of
% material with inner radius degInner and outer radius degOuter, and adding
% those radii to the overall rotor dimensions (rim).

global rotor

fpos = find((rim > Floc),1); % index location of for the outer radius of the failed rim
innerPos = fpos - 1; % position of the failed rim in the rotor

ringInner = Floc - dThicc/2; % Inner radius of degraded ring
ringOuter = Floc + dThicc/2; % Outer radius of degraded ring

% Dimension and position of each new rim made from segmenting the origional
innerPos = failRimPos;
ringPos = failRimPos + 1;
outerPos = failRimPos + 2;

innerRadii = [rotor.radii{failRimPos}(1), ringInner];
ringRadii = [ringInner, ringOuter]; % Degrated ring radii. This vector will always be size [1,2]
outerRadii = [ringOuter, rotor.radii{failRimPos}(2)];

rotor.radii(outerPos:end+2) = rotor.radii(failRimPos:end);
rotor.radii{innerPos} = innerRadii;
rotor.radii{ringPos} = ringRadii;
rotor.radii{outerPos} = outerRadii;

delta(outerPos:end+2) = delta(innerPos:end);
delta(innerPos:outerPos) = 0;

% Lable intact and failed rims
rotor.intact = ones(length(rotor.radii));
rotor.intact(ringPos) = 0

rotor.matInd(outerPos) = rotor.matInd(innerPos);
rotor.matInd(ringPos) = rotor.matInd(outerPos);
