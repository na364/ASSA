function [ki,Ei,T] = beamprops(param,value,mass)
%function [ki,Ei,T] = beamprops(param,value,mass)
%
%function to return the properties of the helium beam, returning
%the values in 1/Angstrom, meV and K respectively.
%
%param - specifies the reference from which the properties are calculated
%value - specifices the value of the above parameter
%mass  - is the atomic mass of helium (3 or 4)
%
%param = 'temp', value is in K
%param = 'energy', value is in meV
%param = 'wavevector', value is in 1/Angstrom
%param = 'chess', value is in mV (Eurotherm reading)

% load in the basic parameters and select the correct atomic mass
load_chess_parameters;
if mass==3
    m = SE_amu * 3.01603;
elseif mass==4
    m = SE_amu * 4.00260;
end
        
if strcmp(param,'temp')
    T=value;
    
    Ei_SI=5/2*SE_kB*T;
    Ei=Ei_SI/SE_e*1000;
    
    ki_SI = (2*m*Ei_SI).^0.5/SE_hbar;
    ki = ki_SI/1e10;
end

if strcmp(param,'energy')
    Ei = value;
    Ei_SI = Ei/1000*SE_e;

    ki_SI = (2*m*Ei_SI).^0.5/SE_hbar;
    ki = ki_SI/1e10;

    T = 2/5*Ei_SI/SE_kB;
end

if strcmp(param,'wavevector')
    ki=value;
    ki_SI = value*1e10;
    
    Ei_SI=SE_hbar^2*ki_SI.^2/(2*m);
    Ei=Ei_SI/SE_e*1000;
    
    T = 2/5*Ei_SI/SE_kB;
end

if strcmp(param,'chess')
    disp('Not yet implemented')
    Ei=0;
    ki=0;
    T=0;
end
