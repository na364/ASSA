% This file contains a list of parameters describing the CHeSS
% apparatus
%
% Note: the parameters are all prefixed SE, to avoid mixing
% with any other variables
%
% Update: 10 Feb 2004, 20 Jan 2005

% General parameters (SI)
SE_h = 6.62608e-34;
SE_hbar = 1.05457e-34; % m^2 * kg * s^(-1)
SE_e = 1.60218e-19;
SE_kB = 1.380658e-23; % m^2 * kg * s^(-2) * K^(-1)
SE_amu = 1.66054e-27;
SE_3hemass = SE_amu * 3.01603;  % mass of 3He
SE_gamma = 2.037895e8;          % gyromagnetic ratio for 3He

% General parameters
SE_hbar_meVps = 6.582119514e-13;

SE_NA = 6.022140857e23;

% solenoid field integrals
SE_Beff1 = 7.26e-3;             % field integral in Tm/A
SE_Beff2 = SE_Beff1;            % placeholder, in case we improve calibration
SE_Beff = (SE_Beff1+SE_Beff2)/2;
