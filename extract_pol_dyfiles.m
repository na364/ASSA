function [base_current,alpha1,real_sig,imag_sig,reps,E0, maxI] = extract_pol_dyfiles(filenameWithPath)
% load parameters required for tilted measurements from dy-files
% 

load (filenameWithPath);

base_current=meas.ibase';
real_sig=(meas.mean.Preal);
imag_sig=(meas.mean.Pimag);

% figure
% hold on
% plot(base_current,real_sig,'b')
% plot(base_current,imag_sig,'r')


% Flipping real_imag signals as function of current, to correct for -1*Ivec
% problem:
% real_sig=fliplr(real_sig)';
% imag_sig=fliplr(imag_sig)';
real_sig=real_sig';
imag_sig=imag_sig';

% figure
% hold on
% plot(base_current,real_sig,'b')
% plot(base_current,imag_sig,'r')

reps=meas.numloops;
alpha1=meas.tilt *pi/180;
E0=meas.beam.E0;

maxI= max(meas.mean.A);
