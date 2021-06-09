function [dK,dkz] = k_xfer(ki,gamma,theta)
%function [dK dkz] = k_xfer(ki,gamma,theta)
%
%function to calculate the total momentum transfer on scattering,
%given the primary incident angle and the total scattering angle
%
%ki - beam wavevector
%gamma - primary incident angle (degrees)
%theta - total scattering angle (degrees)

theta=theta*pi/180;
gamma=gamma*pi/180;

dK = ki*(sin(theta-gamma)-sin(gamma));
dkz = ki*(cos(theta-gamma)+cos(gamma));