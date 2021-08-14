function [base_current,alpha1,real_sig,imag_sig,reps,E0] = extract_pol_dyfiles(meas)


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
if meas.tilt==315&&(strcmp(meas.psu,'deltafine')||strcmp(meas.psu,'deltacoarse'))
    alpha1=135*pi/180; 
end
% when the psu is set to 'deltafine' or 'deltacoarse', the tilt angle is
% set to 315 degrees. However, in our script it should be treated as 135.

E0=meas.beam.E0;

end
