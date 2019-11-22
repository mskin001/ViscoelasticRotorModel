function [newRotor, degMatProp] = degradeRotor(Floc, dThicc)
% This function takes in the failure location and adds in a ring of degraded
% material with a thickness of dThicc and degraded material properties
% degMatProp. This is accomplished by identifying the rim where the failure
% occurs, dividing that into two smaller rims separated by a degraded ring of
% material with inner radius degInner and outer radius degOuter, and adding
% those radii to the overall rotor dimensions (rim).

global mat rim

fInd = find(((rim - Floc) > 1),'first'); % index location of for the inner radius of the failed rim

degInner = Floc - dThicc; % Inner radius of degraded ring
degmOuter = Floc + dThicc; % Outer radius of degraded ring

degRing = [degInner, degOuter] % Degrated ring radii. This vector will always be size [1,2]

newRotor = zeros(1,(length(rim) + length(degRing)));
newRotor(1:fInd) = rim(1:fInd); % All the rim dimentions less than the failure location
newRotor(fInd+1:fInd+2) = degRing; % Add in an additional rim representing the radial failure
newRotor(fInd+3:end) = rim(fInd+1:end); % All rim radii greater than the failure location
