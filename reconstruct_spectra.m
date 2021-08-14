function [lambda_pos,energy,spectrum,corrected_spectrum,exp_lambda,bsf]=reconstruct_spectra(base_current,real_sig,imag_sig,E0,alpha1,tail_current,force_zero, ofac,useWindow)
% This function assumes that current_base is sorted
% alpha1 - tilt angle in radians; alph1=atan(I2/I1)


load_chess_parameters;
m = SE_3hemass;

% manipulate signal
[real_sig,imag_sig] = manipulate_signal(base_current,real_sig,imag_sig,tail_current,useWindow);

real_sig = real_sig(:)';
imag_sig = imag_sig(:)';
base_current = base_current(:)';

% figure(6)
% hold on
% plot(base_current, real_sig, 'b')
% plot(base_current, imag_sig, 'r')

%% Process Time Domain Signal
bs=real_sig+1i*imag_sig;

% if ibase is only positive or negative, mirror
if all(base_current(:) >= 0) || all(base_current(:) <= 0) % is base_current only positive or negative?
    if base_current(2) < 0 % if negative (assuming element #2 would represent that)
        base_current = [base_current fliplr(-base_current(1:end-1))];
        bsm = [bs fliplr(conj(bs(1:end-1)))];
    else
        base_current = [fliplr(-base_current) base_current(2:end)];
        bsm = [fliplr(conj(bs(2:end))) bs];
    end
else
    bsm=bs;
end
N=length(bsm);
% figure; hold on; plot(base_current,real(bsm), 'b'); plot(base_current,imag(bsm), 'r');
%% Wavelength domain
if force_zero==1
    bsm(N)=real(bsm(N));
end
bsf=abs(fftshift(fft(bsm,ofac*N)))*abs(base_current(2)-base_current(1)); % ofac, over sampling
wavelength_res=1/(2*max(base_current)*((SE_gamma/2/pi)*SE_Beff1*m/SE_h))/ofac;
exp_lambda=linspace(-(N*ofac-1)*wavelength_res/2,(N*ofac-1)*wavelength_res/2,N*ofac);

if alpha1==0 % alph1=arctag(I2/I1), when equals zero, we measure a beam profile
    
    exp_lambda=exp_lambda*-1; % the (-1) is here for fixing a specific issue in Cambridge's PSUs (which are connected backwards).
    
    lambda_shift = exp_lambda;
    
else % non beamprofile measurement
    
    K0=beamprops('energy',E0,3)*1e10;
    V0=SE_hbar*K0/m;
    lambda0=SE_h/(m*V0);
    lambda_shift=-lambda0/tan(alpha1)+exp_lambda/sin(alpha1);
    
    if alpha1==pi
        lambda_shift=fliplr(exp_lambda);
    end
    
end


lambda_pos=lambda_shift(lambda_shift>0);
spectrum=bsf(lambda_shift>0);
% figure; 
% plot(exp_lambda(lambda_shift>0),spectrum);
%% Energy domain
energy=SE_h^2./(2*m*lambda_pos.^2);% in joules
jacobian=SE_h/m^0.5./energy.^1.5;
corrected_spectrum=spectrum.*jacobian;
energy=6.2415e21*energy;% convert to meV

end


function [real_sig,imag_sig] = manipulate_signal(base_current,real_sig,imag_sig,tail_current,useWindow)
% Define which section of the ISF is to be used as a background
[~,tail_index]=min(abs(base_current-tail_current));
background_imag=mean(imag_sig(tail_index:end));
imag_sig=imag_sig-background_imag;
real_sig=real_sig-background_imag;

if useWindow
    window=hann(length(base_current));
    real_sig=real_sig.*window;
    imag_sig=imag_sig.*window;
end

end
