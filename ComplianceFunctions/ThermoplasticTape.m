function [S22] = ThermoplasticTape()
% THis function produces the time dependent compliance of the fiberglass
% thermoplastic tape used to make the composite pipes for shockcore. The
% raw data was made availabe by Ha et al. during his masters in 2019.
% Hoop strain vs time was collected during the experiments. The hoop stress
% was calculated using the applied internal pressure, and the thin walled
% pressure vessle assumption. 