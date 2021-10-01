function [base_current,alpha1,real_sig,imag_sig,reps,E0] = extract_pol_dyfiles(meas,intpI)

if exist('intpI','var')&&~isempty(intpI)
    base_current=intpI;
    real_sig=interp1(meas.ibase,meas.mean.Preal,intpI);
    imag_sig=interp1(meas.ibase,meas.mean.Pimag,intpI);
% when the current is not equally spaced, use linear interpolation to make
% the base current equally spaced, so that fft can be used later in reconstruct_spectra.m.
else
    base_current=meas.ibase;
    real_sig=(meas.mean.Preal);
    imag_sig=(meas.mean.Pimag);
end

% figure
% hold on
% plot(base_current,real_sig,'b')
% plot(base_current,imag_sig,'r')

if size(base_current,1)<size(base_current,2)
    base_current=base_current';
end
if size(real_sig,1)<size(real_sig,2)
    real_sig=real_sig';
end
if size(imag_sig,1)<size(imag_sig,2)
    imag_sig=imag_sig';
end

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
