% This scripts aim to calculate the SKW with the same number of points than
% P_lambda, i.e. no interpolation

% load physical constants:
 load_chess_parameters;


%% go from polarization to reconstructed wavelength distribution.
[P_lambda,lambda_m,~,~,~,~] = Full_loop_func.Pkappa2Plambda_single_theta(kappa_vect_s,alpha1,lambda0_m,res_pol_s, Pol_complex_s,list_s,meas.CS_s.Method);

%% Calculate the gridding in energy [J]:

 % initial wavelength
lambda0_m =  1.8e-10;
 
% calculating the range of energies corresponding to the range of wavelengths:
deltaE_J = SE_h^2/(2*SE_3hemass)*(1./lambda_m.^2-1/(lambda0_m)^2);


%% Calculate the SKW from P_lambda and lambda_m

% The Jacobian to go from Plambda to SKw, in [m/J] 
Jacob_det = ...
    SE_h*lambda0_m^3*SE_3hemass./(2*SE_3hemass*deltaE_J*(lambda0_m)^2+(SE_h^2)).^(1.5);


% S(theta_i,DeltaE). Note that this is not properly the
% scattering function because it depends on theta_i and not on
% the momentum transfer! The units are [1/m.J.rad^2]
SKw_tot = P_lambda.*Jacob_det;

