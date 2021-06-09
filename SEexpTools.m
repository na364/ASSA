    
%class spin echo experiment tools: plots and projects WIM
classdef SEexpTools
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        % there are no properties
    end
    
    
    
    
    % define the functions, one by one
    % This class contains
    % - calcMeasurementsParams
    % - projectByTilt
    % - OneD2finalLambda
    % - wavelength2energy
    % - calcMinWithForModeInTilt
    % - modeModelHandle
    % - exploreMeasurements
    % - plotModels
    methods(Static)
        
        
        
        %---------------------------------------------------------------
        % define calcMeasurementsParams function
        %---------------------------------------------------------------
        function [lambda_i_Mat, lambda_f_Mat,wavelengthIntMat, tilt] = calcMeasurementsParams(E0, FWHM, max_dEexchange, theta_i, theta_tot)
            %calcMeasurementsParams to calculate the measurements needed to
            %map a specific mode.
            % modeModelHandle - handle to function which gets a dK, and gives dE.
            % E0 - the beams mean initial energy (in meV)
            
            % what does this do
            load_chess_parameters;
            
            % define the number of points to use 
            numOfPoints = 1e3;
            
            %--------------------------------------------------------------
            %% find beam params (energy, wave-vector, profile)
            %--------------------------------------------------------------
            % incident wavelength matrix
            lambda_i_tmp = linspace(energy2wavelength(E0+2*FWHM,3), energy2wavelength(E0-2*FWHM,3), numOfPoints);
            
            
            
            % calculate the energy of a helium-3 particle with wavelength
            % lambda_i_tmp
            Ei_vec = SEexpTools.wavelength2energy(lambda_i_tmp,3);    
            
            %calculate the wavevector of a helium-3 particle with energy
            %Ei_vec
            ki_vec = energy2wavevector(Ei_vec,3);
            
            
            % convert FWHM of particle energies to standard deviation
            c_param = FWHM/(2*sqrt(2*log(2)));
            
            % something gaussian
            beamIntensity = exp(-(E0-Ei_vec).^2./(2*c_param.^2));
            
            
            % scattered(final, f) wavelength matrix
            lambda_f_tmp = linspace(energy2wavelength(max(Ei_vec)+max_dEexchange,3),energy2wavelength(max(min(Ei_vec)-max_dEexchange,1e-1),3),numOfPoints);            
            
            
            Ef_vec = SEexpTools.wavelength2energy(lambda_f_tmp, 3);
            kf_vec = energy2wavevector(Ef_vec,3);
            
            [lambda_i_Mat, lambda_f_Mat]= meshgrid(lambda_i_tmp,lambda_f_tmp);
            [ki_Mat, kf_Mat] = meshgrid(ki_vec,kf_vec);
            [Ei_Mat, Ef_Mat] = meshgrid(Ei_vec,Ef_vec);
            
            
            
            theta_i_Mat = theta_i + (rand(size(Ei_Mat))-0.5)*0.2;
            theta_out_Mat = theta_tot - theta_i_Mat;
            
            
            
            dE_Mat = Ef_Mat-Ei_Mat;
            dK_Mat = kf_Mat.*sin(theta_out_Mat/180*pi)-ki_Mat.*sin(theta_i_Mat/180*pi);
            %--------------------------------------------------------------
            %--------------------------------------------------------------
            
            
            %--------------------------------------------------------------
            %% find (dk, de) points of phonon
            %--------------------------------------------------------------
            modeModelHandle = SEexpTools.phononModel('SpecBroad');
            modeDE_SpecBroad = modeModelHandle(dK_Mat);
            modeModelHandle = SEexpTools.phononModel('RW');
            modeDE_RW = modeModelHandle(dK_Mat);
      %       modeModelHandle = SEexpTools.phononModel('WaterOnNi6meV');
      %       modeDE_W6 = modeModelHandle(dK_Mat);
%            modeModelHandle = SEexpTools.phononModel('WaterOnNi2meV');
%            modeDE_W2 = modeModelHandle(dK_Mat);
%            modeModelHandle = SEexpTools.phononModel('NiMagnons');
%            modeDE_NiMagnons = modeModelHandle(dK_Mat);
%             modeModelHandle = SEexpTools.phononModel('LongitudinalWaveNi111');
%             modeDE_LA = modeModelHandle(dK_Mat);
%             modeModelHandle = SEexpTools.phononModel('AcousticSurfacePlasmons');
%             modeDE_ASP = modeModelHandle(dK_Mat); %modeDE_ASP(find(abs(modeDE_ASP) > 15))=0;
            
            %find the points in which the scan line crosses the modes (modeDE and -modeDE)
            energyMatchMat = zeros(size(dE_Mat));            
            energyMatchMat = energyMatchMat | abs(modeDE_SpecBroad-dE_Mat) < 1e-1 | abs(modeDE_SpecBroad+dE_Mat) < 1e-1;
            energyMatchMat = energyMatchMat | abs(modeDE_RW-dE_Mat) < 1e-1 | abs(modeDE_RW+dE_Mat) < 1e-1;
     %        energyMatchMat = energyMatchMat | abs(modeDE_W6-dE_Mat) < 1e-1 | abs(modeDE_W6+dE_Mat) < 1e-1;
             %energyMatchMat = energyMatchMat | abs(modeDE_NiMagnons-dE_Mat) < 1e-1 | abs(modeDE_NiMagnons+dE_Mat) < 1e-1;
%           energyMatchMat = energyMatchMat | abs(modeDE_LA-dE_Mat) < 1e-1 | abs(modeDE_LA+dE_Mat) < 1e-1;
%             energyMatchMat = energyMatchMat | abs(modeDE_W2-dE_Mat) < 1e-1 | abs(modeDE_W2+dE_Mat) < 1e-1;
            %energyMatchMat = energyMatchMat | abs(modeDE_ASP-dE_Mat) < 1e-1 | abs(modeDE_ASP+dE_Mat) < 1e-1;

            beamIntensityMat = repmat(beamIntensity,size(energyMatchMat,1),1);
            wavelengthIntMat = energyMatchMat .* beamIntensityMat;
            
            figure; hold on            
            plot(reshape(dK_Mat,1,[]),reshape(modeDE_SpecBroad,1,[]),'b.',reshape(dK_Mat,1,[]),-reshape(modeDE_SpecBroad,1,[]),'b.','MarkerSize',2)
      %      plot(reshape(dK_Mat,1,[]),reshape(modeDE_RW,1,[]),'b.',reshape(dK_Mat,1,[]),-reshape(modeDE_RW,1,[]),'b.','MarkerSize',2)
%             plot(reshape(dK_Mat,1,[]),reshape(modeDE_LA,1,[]),'b.',reshape(dK_Mat,1,[]),-reshape(modeDE_LA,1,[]),'b.','MarkerSize',2)
      %      plot(reshape(dK_Mat,1,[]),reshape(modeDE_W6,1,[]),'b.',reshape(dK_Mat,1,[]),-reshape(modeDE_W6,1,[]),'b.','MarkerSize',2)
%              plot(reshape(dK_Mat,1,[]),reshape(modeDE_W2,1,[]),'b.',reshape(dK_Mat,1,[]),-reshape(modeDE_W2,1,[]),'b.','MarkerSize',2)
%             plot(reshape(dK_Mat,1,[]),reshape(modeDE_NiMagnons,1,[]),'b.',reshape(dK_Mat,1,[]),-reshape(modeDE_NiMagnons,1,[]),'b.','MarkerSize',2)
            %plot(reshape(dK_Mat,1,[]),reshape(modeDE_ASP,1,[]),'b.',reshape(dK_Mat,1,[]),-reshape(modeDE_ASP,1,[]),'b.','MarkerSize',2)
            plot(reshape(dK_Mat,1,[]),reshape(dE_Mat,1,[]),'g.','MarkerSize',2)
            xlabel("$\Delta K$",'Interpreter','latex')
            ylabel("$\Delta E$",'Interpreter','latex')
            
            
            
            figure; hold on
            pcolor(lambda_i_Mat, lambda_f_Mat,wavelengthIntMat); shading flat
            
            
            tmpIndx=find(wavelengthIntMat~=0);
            maxNonZeroLambdaFinal = max(lambda_f_Mat(tmpIndx));
            minNonZeroLambdaFinal = min(lambda_f_Mat(tmpIndx));
            maxNonZeroLambdaIn = max(lambda_i_Mat(tmpIndx));
            minNonZeroLambdaIn = min(lambda_i_Mat(tmpIndx));
            axis([minNonZeroLambdaIn*0.9 maxNonZeroLambdaIn*1.1 minNonZeroLambdaFinal*0.9 maxNonZeroLambdaFinal*1.1])
            
            featureNum = inputdlg('how many features there are?'); featureNum = str2num(featureNum{:});
            zvl = questdlg('Choose a point above the features, then points which seperate the different features in lambda_out, then a point below the lowest feature');
            [l1,l2] = ginput(featureNum+1);
            
            for i=2:length(l2)
                indx = find(lambda_f_Mat < l2(i-1) & lambda_f_Mat > l2(i));
                l1current = lambda_i_Mat(indx);
                l2current = lambda_f_Mat(indx);
                
                intensityCurrent = wavelengthIntMat(indx);        
                indx = find(intensityCurrent > 0.05);                
                p = polyfit(l1current(indx),l2current(indx),1);
                hold on; plot([min(l1current) max(l1current)],p(1)*[min(l1current) max(l1current)]+p(2),'g')
                tilt(i-1) = atan(p(1))*180/pi + 90;
            xlabel("$\lambda_i$")
            ylabel("$\lambda_f$")
            end
            %--------------------------------------------------------------
            %--------------------------------------------------------------
            
        end
        %---------------------------------------------------------------
        %---------------------------------------------------------------
        
        
        
      
        
        
        
        
        %---------------------------------------------------------------
        % define projectByTilt function
        %---------------------------------------------------------------
        function [lambda_1D, lambda_axis_shifted, wavelengthInt_proj,energy,spectrum_in_energy] = projectByTilt(lambda_i_Mat, lambda_f_Mat, wavelengthIntMat_orig,tilt, E0)
            
            lambda_axis_shifted=zeros(size(lambda_i_Mat,2),length(tilt));
            energy=zeros(size(lambda_i_Mat,2),length(tilt));
            spectrum_in_energy=zeros(size(lambda_i_Mat,2),length(tilt));
            
            if ~exist('SE_h','var'), load_chess_parameters; end
            h=figure;
            for i=1:length(tilt)
            
                R=[cosd(tilt(i)-90) -sind(tilt(i)-90) ; sind(tilt(i)-90) cosd(tilt(i)-90)];
                lambda_i_Mat_rotated = R(1)*lambda_i_Mat + R(2)*lambda_f_Mat;
                lambda_f_Mat_rotated = R(3)*lambda_i_Mat + R(4)*lambda_f_Mat;
                F = TriScatteredInterp(lambda_i_Mat_rotated(:),lambda_f_Mat_rotated(:),wavelengthIntMat_orig(:));
                indx = find(wavelengthIntMat_orig > 0.05);
                minMax = [min(min(lambda_i_Mat_rotated(indx)))-1 max(max(lambda_i_Mat_rotated(indx)))+1 min(min(lambda_f_Mat_rotated(indx)))-1 max(max(lambda_f_Mat_rotated(indx)))+1];
                [lambda_i_Mat_new, lambda_f_Mat_new] = meshgrid(linspace(minMax(1),minMax(2),size(lambda_i_Mat,1)),linspace(minMax(3),minMax(4),size(lambda_i_Mat,2)));  
                wavelengthIntMat_rotated = F(lambda_i_Mat_new,lambda_f_Mat_new);

                indx = isnan(wavelengthIntMat_rotated);
                wavelengthIntMat_rotated(indx) = 0;
                wavelengthInt_proj = sum(wavelengthIntMat_rotated,2);
                lambda_1D = lambda_f_Mat_new(:,1);
                
                
                figure; plot(lambda_1D,wavelengthInt_proj); title(['tilted projection peasurement for ' num2str(tilt(i))])
                xlabel("projected wavelength")
                ylabel("projected intensity")
                
                
                [lambda_axis_shifted_tmp, energy_tmp, spectrum_in_energy_tmp] = SEexpTools.OneD2finalLambda(lambda_1D, wavelengthInt_proj, E0, tilt(i));
                lambda_axis_shifted(i,1:length(lambda_axis_shifted_tmp)) = lambda_axis_shifted_tmp';
                energy(i,1:length(energy_tmp)) = energy_tmp';
                spectrum_in_energy(i,1:length(spectrum_in_energy_tmp)) = spectrum_in_energy_tmp';
                
                figure(h)
                subplot(1+ceil(length(tilt)/2),2,i)
                plot(energy(i,:)-E0,spectrum_in_energy(i,:))
                xlabel('dE [meV]'); ylabel('Intensity'); title(['tilt=' num2str(tilt(i))]);
                axis([-10 10 -1e-22 max(spectrum_in_energy(i,:))*1.2])
            end
            figure(h)
            subplot(1+ceil(length(tilt)/2),2,(1+ceil(length(tilt)/2))*2)
            subplot(1+ceil(length(tilt)/2),2,(1+ceil(length(tilt)/2))*2-1)
            
        end
        %---------------------------------------------------------------
        %---------------------------------------------------------------
        
        
        
        
        
        
        
        
        
        
        
        
        
        %---------------------------------------------------------------
        % define OneD2finalLambda function
        %---------------------------------------------------------------
        function [lambda_axis_shifted, energy, spectrum_in_energy] = OneD2finalLambda(lambda_1D, wavelengthInt_proj, E0, tilt)
            
            %
            if ~exist('SE_h','var'), load_chess_parameters; 
            end
            
            
            % find the energy of a helium 3 particle with energy E0
            lambda0=energy2wavelength(E0, 3);
            
            % tilted projection angle (?)
            t1=tilt-90;
            
            
            lambda_axis_shifted=lambda0*tand(t1)+lambda_1D/cosd(t1);
            
            lambda_pos=lambda_axis_shifted(lambda_axis_shifted>0);
            
            spectrum_in_lambda=wavelengthInt_proj(lambda_axis_shifted>0);
            
            energy=SEexpTools.wavelength2energy(lambda_pos,3);
            
            jacobian=(SE_h/SE_3hemass^2)*(2*energy/SE_3hemass).^(-3/2);
            
            spectrum_in_energy=spectrum_in_lambda.*jacobian;
        end
        %---------------------------------------------------------------
        %---------------------------------------------------------------
        
        
        
        
        
        
        
        %---------------------------------------------------------------
        % define wavelength2energy function
        %---------------------------------------------------------------
        function E = wavelength2energy(lambdaInAnrstrem, mass)
            % calculate the energy E of a particle with mass and wavelength
            
            
            % load in the basic parameters and select the correct atomic mass
            if ~exist('SE_amu','var'), load_chess_parameters; end
            
            % if mass is 3 amu i.e. if particle is helium 3
            % converts mass to kg?
            if mass==3
                m = SE_3hemass;
                
            % if mass is 4 amu i.e. if particle is helium 4
            % converts mass to kg?    
            elseif mass==4
                m = SE_amu * 4.00260;
                
            % other masses aren't allowed
            else
                disp('Only masses 3 and 4 allowed')
                return
                
            end
            
            
            % wavelength of particle in SI units: e-10 due to angstroms
            lambda_SI = lambdaInAnrstrem*1e-10;
            
            % energy of the particle from Planck's equation
            E_SI = (SE_h./lambda_SI).^2/(2*m);
            
            %convert to joules (?)
            E = E_SI/SE_e*1000;
        end
        %---------------------------------------------------------------
        %---------------------------------------------------------------
        
        
        
        
        
        
        %---------------------------------------------------------------
        % define calcMinWithForModeInTilt function
        %---------------------------------------------------------------
        function [dK, dE, modeWidthDue2beam] = calcMinWidthForModeInTilt(modeModelHandle, E0, FWHM, max_dK, tiltDeg)
            dK = -max_dK:0.001:max_dK; % in 1/Angstrem, TODO: implement angle resolution of instrument
            dE = modeModelHandle(dK);
            transMat = [cosd(180-tiltDeg), sind(180-tiltDeg)];
            
            for i=1:length(dK)
                lambda_0 = [energy2wavelength(E0,3), energy2wavelength(E0+dE(i),3)];
                lambda_1 = [energy2wavelength(E0+FWHM/2,3), energy2wavelength(E0+FWHM/2+dE(i),3)];                
                modeWidthDue2beam(i) = lambda_0*transMat' - lambda_1*transMat';
            end
            
            plot3(dK,dE,modeWidthDue2beam,'o')
            
        end
        %---------------------------------------------------------------
        %---------------------------------------------------------------
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        %---------------------------------------------------------------
        % define modeModelHandle function
        %---------------------------------------------------------------
        function modeModelHandle = phononModel(modelName)
            switch modelName
                case 'RayleighWaveNi111'
                    modeModelHandle=@(x) 16.88*sin(pi/2/(2.91/2).*x)-0.5192*sin(pi/2/(2.91/2).*x).^3;
                    %                 case 'LongitudinalWaveNi111'
                    %                     modeModelHandle=@(x) *sin(.*x)-*sin(.*x).^3;
                case 'WaterOnNi6meV'
                    modeModelHandle=@(x) (14.53.*(abs(x)<0.6).*abs(x).^2+5.559).*sign(x);
                case 'specular'
                    modeModelHandle=@(x) zeros(size(x));
                case 'WaterOnNi2meV'
                    modeModelHandle= @(x) repmat(2,size(x));
                case 'NiMagnons'
                    modeModelHandle= @(x) (364.*(abs(x)<0.2).*x.^2).*sign(x);
                    %                 case 'AcousticSurfacePlasmons'
                    %                     modeModelHandle=@(x) 4250*x.*(abs(4250*x)<20)+20*(abs(4250*x)>=20);
                case 'SpecBroad'
                    modeModelHandle= @(x) zeros(size(x));
                case 'RW'
                    modeModelHandle= @(x) 17.4*abs(x);
                otherwise
            end
            
        end
        %---------------------------------------------------------------
        %---------------------------------------------------------------
        
        
        
        
        
        
        
        
        
        
        
        
        %---------------------------------------------------------------
        % define exploreMeasurements function
        %---------------------------------------------------------------
        function exploreMeasurements(path2file, filenames)
            
            for i=1:length(filenames)
                
                if mod(i,8) == 1
                    figure; j=1;
                else
                    j=j+1;
                end
                subplot(4,2,j); hold on
                pathAndFilename = [path2file filenames{i} '.mat'];                
                load(pathAndFilename);
                fileTitle=[filenames{i} ', tilt=' num2str(meas.tilt) ', gamma=ga-' num2str(meas.endStatus.gammaSpecular-meas.endStatus.gamma) ...
                           ',azimuth=<110>, Temperature=' num2str(meas.endStatus.tSample) ', E_{beam}=' num2str(meas.beam.E0) ...
                           ',FWHM=' num2str(meas.beam.FWHM) ', 1ML H_2O/Au(111)'];

                [base_current,alpha1,real_sig,imag_sig,reps,E0]=extract_pol_dyfiles(pathAndFilename);
                [lambda_pos,energy,spectrum,corrected_spectrum]=reconstruct_spectra(base_current,real_sig,1*imag_sig,E0,alpha1*180/pi,1.4,0,1);
                plot(energy-E0,corrected_spectrum,'r')
                %plot(lambda_pos,spectrum,'r')
                [lambda_pos,energy,spectrum,corrected_spectrum]=reconstruct_spectra(base_current,real_sig,0*imag_sig,E0,alpha1*180/pi,1.4,0,1);
                plot(energy-E0,corrected_spectrum,'b')
                %plot(lambda_pos,spectrum,'b')
                [lambda_pos,energy,spectrum,corrected_spectrum]=reconstruct_spectra(base_current,0*real_sig,1*imag_sig,E0,alpha1*180/pi,1.4,0,1);
                plot(energy-E0,corrected_spectrum,'m')
                %plot(lambda_pos,spectrum,'m')
                axis([-10 10 -1e-23 1e-21])
                legend('abs','real','imag')
                title(fileTitle)
                %xlabel('dE [meV]'); ylabel('Intensity')
            end
        end
        %---------------------------------------------------------------
        %---------------------------------------------------------------
        
        
        
        
        
        
        
        
        
        
        
        %---------------------------------------------------------------
        %define plotModels function
        %---------------------------------------------------------------
        function plotModels()
            dK = -10:0.05:10;
            modeModelHandle1 = SEexpTools.phononModel('RayleighWaveAu111');
            modeModelHandle2 = SEexpTools.phononModel('WaterOnAu6meV');
            modeModelHandle3 = SEexpTools.phononModel('WaterOnAu2meV');
            modeModelHandle4 = SEexpTools.phononModel('LongitudinalWaveAu111');
            modeModelHandle5 = SEexpTools.phononModel('AcousticSurfacePlasmons');
            plot(dK,modeModelHandle1(dK),'b',dK,-modeModelHandle1(dK),'b'); hold on
            plot(dK,modeModelHandle2(dK),'b',dK,-modeModelHandle2(dK),'b'); hold on
            plot(dK,modeModelHandle3(dK),'b',dK,-modeModelHandle3(dK),'b'); hold on
            plot(dK,modeModelHandle4(dK),'b',dK,-modeModelHandle4(dK),'b'); hold on
            plot(dK,modeModelHandle5(dK),'b',dK,-modeModelHandle5(dK),'b'); hold on            
        end
        %---------------------------------------------------------------
        %---------------------------------------------------------------
        
        
        
        
        
        
   
    end
    
end

