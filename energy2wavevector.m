

function [k] = energy2wavevector(E, mass)
%function [k] = energy2wavevector(E, mass)
%
%function to convert the energy of atoms in a helium beam to the
%associated wavevector
%
%k is in 1/Angstrom
%E is in meV
%mass is either 3 or 4




% if the parameters exist, use them
if ~exist('SE_amu','var'), load_chess_parameters; 

%end condition
end



% if mass is 3 amu i.e. if particle is helium 3, then find mass in SI
if mass==3
    m = SE_amu * 3.01603;
    
% if mass is 4 amu i.e. if particle is helium 4, then find mass in SI
elseif mass==4
    m = SE_amu * 4.00260;
    
end



%find energy in SI
Ei_SI = E/1000*SE_e;

% find wavevector in SI
% in inverse angstroms
ki_SI = (2*m*Ei_SI).^0.5/SE_hbar;

%convert to inverse metres
k = ki_SI/1e10;
