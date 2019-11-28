function [delta] = degradeRotor(rim, Floc, dThicc, degStiffPerc, mat, delta)
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

% Parse failed rim into damaged and undamaged sections - treated as new rims
rotor.pos(failRimPos:end) = rotor.pos(failRimPos:end)+2;
rotor.pos(end+1:end+2) = [failRimPos, ringPos];
rotor.pos = sort(rotor.pos);

rotor.radii(newOuterPos:end+2) = rotor.radii(failRimPos:end);
rotor.radii{newInnerPos} = newInnerRadii;
rotor.radii{ringPos} = ringRadii;
rotor.radii{newOuterPos} = newOuterRadii;

rotor.Q{newOuterPos} = rotor.Q{newInnerPos};
rotor.rho(newOuterPos) = rotor.rho(newInnerPos);
rotor.stren{newOuterPos} = rotor.stren{newInnerPos};

delta(newOuterPos:end+2) = delta(newInnerPos:end);
delta(newInnerPos:newOuterPos) = 0;

% Assign material properties for the degraded ring
baseProps = load(mat.file{failRimPos});
baseProps.mstiff(2) = degStiffPerc * baseProps.mstiff(2);
rotor.Q{ringPos} = stiffMat(baseProps.mstiff,'no'); % 'no' to disable VE part of stiffMat
rotor.rho(ringPos) = baseProps.rho;
rotor.stren{ringPos} = baseProps.stren;
