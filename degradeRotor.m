function [newRotor, rimInd] = degradeRotor(Floc, dThicc)
% This function takes in the failure location and adds in a ring of degraded
% material with a thickness of dThicc and degraded material properties
% degMatProp. This is accomplished by identifying the rim where the failure
% occurs, dividing that into two smaller rims separated by a degraded ring of
% material with inner radius degInner and outer radius degOuter, and adding
% those radii to the overall rotor dimensions (rim).

global mat rim

rimInd = find(((rim - Floc) > 1),'first'); % index location of for the inner radius of the failed rim

degInner = Floc - dThicc/2; % Inner radius of degraded ring
degOuter = Floc + dThicc/2; % Outer radius of degraded ring

degRing = [degInner, degOuter]; % Degrated ring radii. This vector will always be size [1,2]

newRotor = zeros(1,(length(rim) + length(degRing)));
newRotor = [rim(1:rimInd), degRing, rim(rimInd+1:end)]; % All the rim dimentions less than the failure location
