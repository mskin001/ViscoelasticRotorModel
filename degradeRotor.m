function [newRotor, rimInd] = degradeRotor(rim, Floc, dThicc)
% This function takes in the failure location and adds in a ring of degraded
% material with a thickness of dThicc and degraded material properties
% degMatProp. This is accomplished by identifying the rim where the failure
% occurs, dividing that into two smaller rims separated by a degraded ring of
% material with inner radius degInner and outer radius degOuter, and adding
% those radii to the overall rotor dimensions (rim).

global rotor

fpos = find(((rim - Floc) > 0),1); % index location of for the outer radius of the failed rim
failRimPos = fpos - 1; % position of the failed rim in the rotor

ringInner = Floc - dThicc/2; % Inner radius of degraded ring
ringOuter = Floc + dThicc/2; % Outer radius of degraded ring

% Dimension and position of each new rim made from dividing the origional failed rim.
newInnerRadii = [rotor.radii{failRimPos}(1), ringInner];
newInnerPos = failRimPos;
ringRadii = [ringInner, ringOuter]; % Degrated ring radii. This vector will always be size [1,2]
ringPos = failRimPos + 1;
newOuterRadii = [ringOuter, rotor.radii{failRimPos}(2)];
newOuterPos = failRimPos + 2;

rotor.pos(failRimPos:end) = rotor.pos(failRimPos:end)+2;
rotor.pos(end+1:end+2) = [failRimPos, ringPos];
[rotor.pos,ind] = sort(rotor.pos);


rotor.Q{ringPos+1} = rotor.Q{ringPos-1};
rotor.rho(ringPos+1) = rotor.rho(ringPos-1);
rotor.stren{ringPos+1} = rotor.stren{ringPos-1};