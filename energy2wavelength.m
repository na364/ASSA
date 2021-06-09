% function lambda = energy2wavelength(energy,mass)
%
% function to convert energy (in meV) to wavelength (in Angstrom)
% wavelength is returned in Angstroms




function lambda = energy2wavelength(energy, mass)

% h is Planck's constant
h = 6.62608e-34;

% e is the charge on the electron
e = 1.60218e-19;

%c is the speed of light
c = 2.99792e8;


% load in the basic parameters and select the correct atomic mass
if ~exist('SE_amu','var'), load_chess_parameters; 

% end condition    
end




% if mass is 3 amu i.e. if particle is helium 3
% converts mass to kg?
if mass==3
    m = SE_amu * 3.01603;


% if mass is 4 amu i.e. if particle is helium 4
% converts mass to kg?
elseif mass==4
    m = SE_amu * 4.00260;

    
% helium only has two isotopes - anything else is not allowed
else
    disp('Only masses 3 and 4 allowed')
    
    %return nothing
    return
    
% end conditions    
end




% calculate the energy in SI units
% change from eV to SI?
E_SI = energy ./ 1000*SE_e;

% calculate wavelength in metres
% this is the output
% multiply by e10 to cancel out the angstrom-based input
lambda = SE_h ./ sqrt(2*m*E_SI) * 1e10;

