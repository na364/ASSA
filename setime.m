function [result] = setime(param,value,ki)
%function [result] = setime(param,value,ki)
%
%function to convert between spinecho time and solenoid current,
%given the calibration for the coils
%
%param="current", value is in amps, result is in ps
%param="time", value is in ps, result is in amps
%ki is the beam wavevector

load_chess_parameters
ki_SI = ki*1e10;
a = 4*pi^2*SE_gamma*SE_3hemass^2/(SE_h^2*ki_SI.^3)*SE_Beff1;

if strcmp(param,'current')  
    tse_SI = a*value;
    result = tse_SI/1e-12;
end
    
if strcmp(param,'time')  
    tse_SI = value*1e-12;
    result = tse_SI/a;
end
