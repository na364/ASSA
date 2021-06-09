classdef Full_loop_func
    % FULL_LOOP_FUNC is the library of functions which allow to perform
    % the full loop from the polarization to the intermediate scattering
    % function in the compressed sensing method
    
     properties(Constant)
     end
     
     methods(Static)
         
         %% The main function:
         function [meas] = CS_scan_current(varargin)
                % CS_SCAN_CURRENT is the function which calls scan_current
                % and Pol2Isf where the CS is performed 
                %
                % The arguments are:
                %       psu:                                               select 'deltafine', 'deltacoarse'
                %
                %       ibase_matrix:                                      matrix with the vector of currents to
                %                                                          measure for each loop
                %
                %       numloops:                                          number of loops
                %
                %       savemode:                                          1 (on) or 0 (off)
                %
                %       CS_s:                                              structure containing the information to  perform compressed sensing. 
                %                                                          It is provided by the function Distrib_func.gen_ibase_matrix
                %
                %            CS_s.ibasetoreconstruct                       vector of currents to reconstruct. Its number of points is the "res" or resolution
                %
                %            CS_s.ibase_matrix                             matrix where each row is the vector of currents to measure in the loop
                %
                %            CS_s.list_mask_matrix                         matrix where each row contains the list of res 1s and 0s corresponding 
                %                                                          to the measured/non measured currents  
                %
                %            CS_s.list_index_matrix                        matrix where each row contains the list of the size(CS_s.list_mask_matrix,1) 
                %                                                          indexes corresponding to the measured currents 
                %
                %            CS_s.METHOD                                   Method for the CS reconstruction. Choose among 'DFT (FT based), 'DWT' (wavelet based),
                %                                                           'Continuous' (continuous method using wavelets)   
                %
                %
                % The ouptut is saved in a file, if savemode is 'on' or
                % 'off'  or in a structure. The structure contains the
                % following fields:
                %      meas:                                               strucutre containing the raw data and the
                %                                                          information about the measurement. See scan_current
                %                                                          for acomplete description
                %      meas.CS_s:                                          structure which contains the CS_s structure as described above
                %      meas.cs_results:                                    structure which contains the CS results and which is delivered by Pol2ISF. A description of the entries
                %                                                          can be found in this function.
                %

                %% READ IN INPUT PARAMETERS AND SET DEFAULTS IF NECESSARY
                for i=1:2:length(varargin)
                    switch varargin{i}
                        case 'psu'
                            psu = varargin{i+1};
                        case 'ibase_matrix'
                            ibase_matrix = varargin{i+1};     
                        case 'numloops'
                            numloops = varargin{i+1};            
                        case 'savemode'
                            savemode = varargin{i+1};                        
                        case 'CS_s'
                            CS_s = varargin{i+1};  
                        case 'file_test'
                            file_test = varargin{i+1};                            
                        otherwise
                            warning(['variable ' varargin{i} ' is not defined']);
                    end
                end

                %% SET DEFAULT PARAMETERS

                % check mandatory variables have been specified
                if exist('psu','var')==0
                    error('You must specify which PSU you are using')
                end

                if exist('numloops','var')==0
                    numloops=1;
                    disp('Using default: numloops=1');
                end

                if exist('ibase_matrix','var')==0
                % ibase defaults differently, for different power supplies.
                % We generate the ibase_matrix and the CS_s structure
                % through Distrib_func.gen_ibase_matrix with minimum
                % parameters and doing a full sampling!!  
                    switch psu
                        case 'deltafine'
                            [ibase_matrix,CS_s] = Distrib_func.gen_ibase_matrix('Imin',-1,'Imax',1,'deltaI',2e-2,'numloops',numloops); 
                            %ibase_matrix=repmat(linspace(-1,1,21),numloops);
                            disp('Using main psu default: Imin = -1 A,Imax = 1 A,deltaI = 1e-2 A');
                        case 'deltacoarse'
                           [ibase_matrix,CS_s] = Distrib_func.gen_ibase_matrix('Imin',-0.01,'Imax',11.3,'deltaI',1e-1,'numloops',numloops); 
                            %ibase_matrix=repmat(linspace(-0.01,11.3,21),numloops);
                            disp('Using main psu default: ibase=linspace(-0.01,11.3,21)');
                        case 'fine'
                            %TODO add default for fine
                    end
                end

                % In the case one wants to provide an ibase_matrix but not
                % a CS_s structure, the CS_s is builded artificially (in
                % particular the list_mask_matrix and list_index_matrix). 
                % However this might not be an optimum reconstruction since
                % ibase_matrix might not be well choosen 
                if exist('CS_s','var')==0
                      if exist('ibase_matrix','var')==1
                           [~,CS_s] = Distrib_func.gen_ibase_matrix('Imin',min(ibase_matrix),'Imax',max(ibase_matrix),'deltaI',1e-3,'numloops',numloops);                      
                           CS_s.ibase_matrix = ibase_matrix;
                           CS_s.list_mask_matrix = zeros(size(ibase_matrix,1),numloops);
                           for kk = 1:numloops
                               CS_s.list_index_matrix(:,kk) = find(CS_sibase2reconstruct == ibase_matrix(:,kk));
                               CS_s.list_mask_matrix( CS_s.list_index_matrix(:,kk),kk) = 1;                           
                           end
                           disp('Generating CS_s from the function Distrib_func.gen_ibase_matrix with standard parameters. But be carefull because maybe your ibase is not well choosen and the reconstruction could not be optimal.');
                            %error('If you want to perform CS you need to provide a CS_s structure.');                          
                      end
                end
                
                 if exist('savemode','var')==0
                    savemode=1;
                    disp('Using default: savemode=1');
                 end
                
                   if exist('file_test','var')==0
                    test_mode=0;
                    disp('No test mode');
                   else 
                     test_mode = 1;
                     disp(['test mode is on on the file'  file_test]);
                  end
                
                 
                 
                %% PREPARE FILE TO SAVE THE DATA

              
                %set up the file
                global data_file;
                
                if (savemode==1)
                        define_datafile_mfile('dy');
                        disp(['Storing data to: ' data_file])
                        load(data_file);
                end


                %% PREPARE FIGURE TO DISPLAY CS RESULTS


                %clear all;
                fh = figure(2);

                if ishandle(fh) % avoids ovewriting figure 2 (= result's GUI)
                    close(fh);
                    fh = figure(2);
                    [fig_s] = Full_loop_func.Create_display_figure(fh);
                end


                % ensure the abort flag is removed
                if exist('/icon/ABORT','file')
                    system('rm /icon/ABORT')
                end
                abort_me = 0;

                %% Set general fields in meas structure
                % not including fields that scan_current calculates (these
                % will be done later on the 1st loop                

                if (test_mode == 0)
                    meas.psu = psu;
                    meas.ibase_matrix = ibase_matrix;                
                    meas.numloops = numloops;                
                    % meas.setime = THIS IS Meaningless;
                    meas.CS_s  = CS_s; % compressed sensing parameters                
                else                
                    meas_test1 = load(file_test);                   
                    meas = meas_test1.meas;
                    meas_orig = meas_test1.meas; % keep original measurement structure for comparison with cs results
                    
                    % simulate undersampling
                    for loop_index = 1:numloops
                        ibaseSampled = CS_s.ibase_matrix(:,loop_index);
                        indx = [];
                        for tmp=1:length(ibaseSampled)
                            indx = [indx find(abs(meas.ibase-ibaseSampled(tmp))<1e-10,1)];
                        end
                        meas.loop(loop_index).Preal = meas_orig.loop(loop_index).Preal(indx);
                        meas.loop(loop_index).Pimag = meas_orig.loop(loop_index).Pimag(indx);
                        meas.loop(loop_index).Pmag = meas_orig.loop(loop_index).Pmag(indx);
                        meas.loop(loop_index).deltaI0 = meas_orig.loop(loop_index).deltaI0(indx);
                        meas.loop(loop_index).deltaPhase = meas_orig.loop(loop_index).deltaPhase(indx);
                    end
                    
                    meas.ibase_matrix = ibase_matrix; 
                    meas.CS_s  = CS_s;
                end
                
                %% LOOP OVER CURRENTS + CS AT THE END OF EACH LOOP:

                for loop_index = 1:numloops

                    % This is the code to pause the loop
                    % Check the measurement has not been paused
                    % if it has, then wait here until the PAUSE file has been removed...

                    if exist('/icon/PAUSE','file')==2
                        disp(['Paused ' data_file ' at the end of loop ' num2str(loop_index)]);
                        while exist('/icon/PAUSE','file')==2
                            pause(1);
                        end
                        disp('Measurement resumed');
                    end

                    % This is the code to abort the loop
                    % check if the measurement has been set to abort at the end
                    % of this loop.  if so, notify on screen and set abort flag
                    if exist('/icon/ABORT','file')==2
                        if abort_me==0
                            disp(['Aborting ' data_file ' at the end of loop ' num2str(loop_index)]);
                            abort_me = 1;
                            meas.numloops = loop_index;
                        end
                    end


                       %% Call scan_current with save mode off
                        if (test_mode == 0)
                        	meas_one_loop = scan_current('psu',psu,'ibase',ibase_matrix(:,loop_index),'numloops',1,'savemode','off');
                        else
                            meas_one_loop = meas;
                        end
                       
                         % Record constants for every loop from loop 1
                        if loop_index == 1
                            meas.tilt = meas_one_loop.tilt;            %tilt should be recorded in degrees
                            meas.count_time = meas_one_loop.count_time;
                            meas.beam = meas_one_loop.beam;
                            meas.phasecoil = meas_one_loop.phasecoil;
                            meas.sepoint = meas_one_loop.sepoint;
                        end

                         % This saves the data of each loop
                         meas.loop(loop_index) = meas_one_loop.loop(loop_index); % raw data

                        % This saves the constants for every loop from last loop
                        if loop_index == meas.numloops
                            meas.endBkg = meas_one_loop.endBkg;
                            meas.endStatus = meas_one_loop.endStatus;                       
                        end

                       
                        
                         % interstitial save
                       if (savemode==1)
                            save(data_file,'meas');
                        end 


                        %% Compressed sensing:

                        % Recombine the data which have been measured up to loop_index 
                        %and provides an updated list of measured current indexes:
                        [meas,list_new] = Full_loop_func.concat_avg_data(meas,loop_index);
                       
                        % Do compressed senssing (from polarization to wavelength distribution) 
                        %and then full loop (up to the ISF) and plot the results:
                        [cs_result_struct] = Full_loop_func.Pol2ISF(meas,list_new);
                        
                         % Plot the CS results:
                       Full_loop_func.update_display_figure('fig_s',fig_s,'cs_result_struct',cs_result_struct,'meas',meas_orig);
                 
                         % Saves the CS results in the last loop
                        if loop_index == meas.numloops
                           meas.cs_results = cs_result_struct;
                        end

                        %% end of loop save
                        if (savemode==1)
                            save(data_file,'meas');
                        end 


                        % abort at end of loop is asked
                        % exit the loop if we've been asked to abort

                        if exist('/icon/ABORT','file')==2
                            if abort_me==0
                                disp(['Aborting ' data_file ' at the end of scan ' num2str(loop_index)]);
                                abort_me = 1;
                                meas.numloops = loop_index; 
                            end
                            system('rm /icon/ABORT');
                            break
                        end


                   
                end


            end
         
         
         
        %% Full loop functions:
        
          % Full loop: from polarization towards the intermediate scattering function:
         function [cs_result_struct] = Pol2ISF(meas,list_new)
            %POL2ISF: This function calculates the loop from the polarization to the
            %"pseudo" intermediate scattering function I(theta,time). The first step consists in applying the compress sensing (CS) method
            % to reconstruct the wavelength distribution from an
            % undersampled polarization measurement
            %
            %
            %   The arguments are the following ones:
            %       MEAS: the structure containing the experimental set-up
            %       and the data. In particular we need the following
            %       entries:
            %           MEAS.MEAN.PREAL: real part of the polarization.
            %           This vector contains all the measured points up to
            %           index_loop
            %           MEAS.MEAN.PIMAG: idem for the imaginary part of
            %           the polarization
            %           MEAS.IBASE_MATRIX: the ibase matrix which contains
            %           the measured currents
            %           MEAS.CS_S.LIST_MASK_MATRIX: the list of 1 and 0s
            %           corresponding to the measured/non-measured currents
            %           in index_loop
            %           MEAS.CS_S.METHOD: the method of reconstruction for
            %           CS, which is one among 'DFT', 'DWT' and
            %           'Continuous'
            %           MEAS.CS_S.IBASETORECONSTRUCT: the ibase vector to
            %           reconstruct
            %         
            %       LIST_NEW: updated list of measured currents
            %
            %   The output of the function are stored in the following structure:
             %      CS_RESULT_STRUCT
            %                CS_RESULT_STRUCT.LAMBDA_M: the wavelength vector in [m], 
            %                CS_RESULT_STRUCT.P_LAMBDA: the distribution of wavelengths in  [1/m^2*rad^2]
            %                CS_RESULT_STRUCT.ENERG_MEV: the uniformely spaced energy
            %                grid in [meV]
            %                CS_RESULT_STRUCT.FREQ_HZ: the uniformely spaced frequency
            %                grid in [Hz]
            %                CS_RESULT_STRUCT.SKW: the "pseudo" scattering function for a
            %                single theta angle in [1/(J*m*rad)]
            %                CS_RESULT_STRUCT.TIME_PS: the uniformely spaced time vector
            %                in [ps]
            %                CS_RESULT_STRUCT.IKT : the "pseudo" intermediate
            %                scattering function in [1/(m)]
            %                CS_RESULT_STRUCT.LAMBDAF_M_INTERP; the
            %                interpolated wavelength vector which is used
            %                to calculate the uniformely spaced energy
            %                gridding. In [m].
            %                CS_RESULT_STRUCT.PLAMBDA_INTERP: the interpolated
            %                wavelength distribution from which we
            %                calculate the SKW (corresponding to the
            %                uniformely spaced energy gridding). In  [1/m^2*rad^2]
            %                CS_RESULT_STRUCT.KAPPA_VECT_INVM: the symmetrical
            %                kappa vector calculated from
            %                meas.mean.CS_s.ibasetoreconstruct
            %                CS_RESULT_STRUCT.PLAMBDA_ORIG: reconstructed
            %                wavelength corresponding to the tilted
            %                wavelength vector. This one is usefull for
            %                retrieving the reconstructed polarization.  In  [1/m^2*rad^2]. 
            %                CS_RESULT_STRUCT.LAMBDAM_ORIG: the wavelength vector directly
            %                calculated from the kappa vector (and from ibasetoreconstrcut). This one is usefull for
            %                normalizing the reconstructed polarization.  In  [1/m]. 
            %                CS_RESULT_STRUCT.POL_RECONSTR: the reconstructed
            %                polarization.
            %                CS_RESULT_STRUCT.CONV: structure with conversion
            %                factors
        

               
                % Load the paths of the CS files:
%                 addpath /icon/scripts/spgl1-1.8;
                  addpath spgl1-1.8;
%                 addpath /icon/scripts/;

                % Conversion constant from Joules to meV:
                C_J2meV = 6.24e21; % in meV/J

                %% Load physical parameters and instrument parameters
                load_chess_parameters;

               

                %% Initialize energy, wavelength and current:      
                 E0_meV = meas.beam.E0;
              
                 [lambda0_m,K0_m,~,conv_factors] = Full_loop_func.calc_E0( E0_meV);


                 % Conversion factor between I and kappa (the Fourier pair of wavelength) (See Eq. 3.3 of Gil's thesis)
                C_I2kappa = SE_gamma*SE_3hemass*SE_Beff/(2*pi*SE_h); %in [1/(A.m)]

                % Tilted kappa vector in [1/(m)]
                kappa_vect =  C_I2kappa*meas.CS_s.ibasetoreconstruct;

                % SE time in [ps]:
                time_vect_ps = setime('current',meas.CS_s.ibasetoreconstruct/sqrt(2),K0_m*1e-10);

                % read complex polarization
                Pol_complex = complex(meas.mean.Preal,meas.mean.Pimag);

                % Mirror polarization and intensity if needed
                %Pol_complex_s = Pol_complex;
                %kappa_vect_s = kappa_vect;
                %ibase_s = analysis.ibase;
                %list_s = list;
                [Pol_complex_s, kappa_vect_s,ibasetoreconstruct_s,list_s] = Full_loop_func.mirror('Pol',Pol_complex,'kappa_vect',kappa_vect,'ibase',meas.CS_s.ibasetoreconstruct,'list',list_new); % mirror the data

                % redefine the "resolution"  ie. the number of points to be
                % reconstructed:
                res_pol_s = length(list_s);



                % tilting angle of current vector in [rad]. Note that if meas.tilt =
                % 315deg we scan the currents from I1_min ->I1_max but if meas.tilt =
                % 135 deg we scan the currents from I1_max ->I1_min. alpha1 should be  
                alpha1 = meas.tilt*pi/180; % in rad

                % Calculate distribution of wavelengths: COMPRESS SENSING
                [P_lambda,lambda_m,Coeff,AdditionalInfo,P_lambda_original,lambda_m_original] = Full_loop_func.Pkappa2Plambda_single_theta(kappa_vect_s,alpha1,lambda0_m,res_pol_s, Pol_complex_s,list_s,meas.CS_s.Method);
               
                % Calculate the reconstructed polarization:
                [Pkappa_reconstr] = Full_loop_func.Plambda2Pkappa_single_theta(P_lambda_original,meas.CS_s.Method,kappa_vect_s);

                % Calculation of a uniformely space gridding of energies
                [DeltaE_grid_unif_J,DeltaE_grid_nonunif_J,lambda_interp_m] = Full_loop_func.lambda2energ(lambda_m,lambda0_m,K0_m,time_vect_ps);
                
                % Calculation of the scattering function on both non-evenly (SKW_tot) and evenly (SKW_tot_interp) spaced grid of (DeltaE,theta_i_grid)_DeltaK.
                [SKw_tot,SKw_tot_interp,~] = Full_loop_func.Plambda2Sthetaw_v2(P_lambda,lambda_interp_m,DeltaE_grid_unif_J,DeltaE_grid_nonunif_J,lambda0_m,alpha1,meas.CS_s.Method,Coeff,AdditionalInfo);

                % The pseudo intermediate scattering function I(t) for a given value of DeltaK and theta_i
                [ IKt_tot,time_vect] = Full_loop_func.Ithetat_calc(SKw_tot_interp,DeltaE_grid_unif_J);
         
                % Copy the results into the output structure
                cs_result_struct.lambda_m            = lambda_m;
                cs_result_struct.Plambda             = P_lambda;
                cs_result_struct.Energ_meV           = DeltaE_grid_nonunif_J*C_J2meV;
                cs_result_struct.Freq_Hz             = (DeltaE_grid_nonunif_J*C_J2meV)*conv_factors.C_J2Hz;
                cs_result_struct.SKw                 = SKw_tot;
                cs_result_struct.time_ps             = time_vect;
                cs_result_struct.IKt                 = IKt_tot;
               % cs_result_struct.lambdaF_m_interp    = lambda_interp_m; 
               % cs_result_struct.Plambda_interp      =  P_lambda_interp_lambda; %it is created in Plambda2Sthetaw but only in the mothod of 'Continuous Wavelet'                
                cs_result_struct.Plambda_orig        = P_lambda_original;
                cs_result_struct.lambdam_orig        = lambda_m_original;
                cs_result_struct.kappa_vect_invm     = kappa_vect_s; 
                cs_result_struct.ibasetoreconstruct  = ibasetoreconstruct_s; 
                cs_result_struct.Pol_reconstr        = Pkappa_reconstr; 
                cs_result_struct.conv = conv_factors;


               

         end

         
         
         % From polarization towards wavelength distribution
         function [P_lambda_cut,lambda2_cut_m,Coeff,AdditionalInfo,P_lambda_orig,lambda_vect_m] = Pkappa2Plambda_single_theta(kappa_vect_rot_m,alpha,lambda0_m,res,P_kappa_slice,list,Method)
            %P_LAMBDA_CALC(KAPPA_VECT_ROT_M,RES,P_KAPPA_SLICE,LIST,FIG_PLOT) Fourier transforms and reconstructs the wavelength
            % distribution along the line of measurement determined by the currents
            % relations, making use of the Compresed Senssing codes of A. Jones. 
            %
            % The arguments are the following ones:
            %   KAPPA_VECT_ROT_M: kappa vector (measurement line), in [1/m]          
            %   ALPHA: direction of the line of measurement : I2/I1 = tan(alpha1) in
            %   [rad]
            %   LAMBDA0_M: incident beam wavelength in [m]
            %   RES: "resolution" or number of points to be reconstructed. This
            %   determines the dimension of the final wavelength distribution. It is a
            %   scalar
            %   P_KAPPA_SLICE: matrix containing the polarization intensity allong the kappa_vect_rot, in
            %   [1/rad^2] for each scattering angle 
            %   LIST: masking vector where 1 corresponds to a measured point
            %   (P_kappa_slice(list(j)==1)!=0) and 0 corresponds to a non-measured
            %   point (_kappa_slice(list(k)==1)=0).
            %   METHOD: method for CS reconstruction. 
            %
            %
            % This function gives back the following outputs:
            %   P_LAMBDA_CUT: matrix containing the final wavelength distribution, in
            %   [1/m^2*rad^2]
            %   LAMBDA2_CUT_M: vector containing the final wavelength (the projection
            %   of the wavelength vector corresponding to the kappa_vector, projected
            %   onto the final wavelength axis). It only contains
            %   the positive (and absolute) values of the final wavelength.  The units are [m].
            %   COEFF: if the CS reconstructed method is Continuous, this
            %   contains the Fourier coefficients (the amplitude of each
            %   component of the Fourier series).
            %   ADDITIONALINFO: additional info for CS
            %   PLAMBDA_ORIG: wavelength distribution obtained from the FT
            %   and reconstruction of the polarization [1/m^2*rad^2]
            %   LAMBDA_VECT_M: wavelength vector orignal obtained from the
            %   current vector. (m)


                %% Wavelengths matrix:

                % Note that because of the  Fourier slice-projection theorem, the 1D FT of
                % P_kappa_onescan gives directely the projected intensity of the wavelength
                % distribution long the axis (lambda_1,lambda_2) perpendicular to the
                % direction defined by (kappa1,kappa2)


                % projection angle (see Figure B.1 of Appendix B of Gil's thesis)
                 alpha_prime = alpha-pi/2;

                % Calculate the gridding in wavelength corresponding to the gridding in
                % kappa.
                kappa_min = min(kappa_vect_rot_m); % in [1/m]
                kappa_max = max(kappa_vect_rot_m); % in [1/m]

                delta_kappa_vect_rot_m = (kappa_max-kappa_min)/length(kappa_vect_rot_m);

                kappa_min_int=round(kappa_min/delta_kappa_vect_rot_m);
                kappa_max_int=round(kappa_max/delta_kappa_vect_rot_m);

                delta_lambda_m = 1/(res*delta_kappa_vect_rot_m); % in [m]
                lambda_min_m = -round(res/2)*delta_lambda_m; % in [m] and relative value (0.1*1e-10-lambda0_m*tan(alpha_prime))*cos(alpha_prime);%

               %============== Compress sensing code: =================================

               % should return lambda_vect in [m] and P_lambda in [1/m^2*rad^2]

               Filter = Compressed_sensing.MakeONFilter('Daubechies',8);
               %Filter=MakeONFilter('Haar',2);
               Levels=5;
               Coeff=1;
               AdditionalInfo=1;
               switch Method
                   case 'DFT'
                    [P_lambda_orig,lambda_vect_m] =...
                        Compressed_sensing.FastFourierReconFromSamplesZeroCentered(P_kappa_slice,list,res,kappa_min,kappa_max);
                   case 'DWT'
                    [P_lambda_orig,lambda_vect_m] =...
                        Compressed_sensing.FastFourierReconFromSamplesZeroCenteredWithWavelets(P_kappa_slice,list,res,kappa_min,kappa_max,Filter,Levels);
                   case 'ContinuousWavelet'
                    [Coeff,AdditionalInfo]=Compressed_sensing.SolveForWaveletCoefficientsFromFourierSamples(P_kappa_slice,list,Filter,[kappa_min_int,kappa_max_int],delta_kappa_vect_rot_m,lambda_min_m);
                        fineres=res;
                        lambda_vect_m=(0:fineres-1)*delta_lambda_m*res/fineres +lambda_min_m;  
                        P_lambda_orig=Compressed_sensing.EvaluateFunctionFromWaveletCoefficients(lambda_vect_m,Coeff,AdditionalInfo);
                   otherwise
                       error('Invalid Reconstruction Method Chosen')
               end


                %Final wavelength, taken from appendix B of gil's thesis
                lambda2_m = lambda_vect_m./cos(alpha_prime)+lambda0_m*tan(alpha_prime); 

                % select only positive wavelengths:
                lambda2_m = lambda2_m(lambda2_m>0);
                P_lambda = P_lambda_orig(lambda2_m>0);

                % Taking a range of positive lambda2 where min(lambda2) is not too small,
                % in order to avoid very large energy transfers
                cutoff = 0; %floor(length(lambda2_m)/20);
                lambda2_cut_m = lambda2_m(cutoff+1:length(lambda2_m)-cutoff); %in [m]

                P_lambda_cut =  P_lambda(cutoff+1:length(P_lambda)-cutoff);
    

         end
         
         % From wavelength distribution towards polarization
         function [Pkappa_reconstr] = Plambda2Pkappa_single_theta(Plambda,Method,kappa_vect_s)
            %PLAMBDA2PKAPPA_SINGLE_THETA Fourier transforms back the reconstructed
            %Plambda into the reconstructed Polarization in kappa units. 
            %
            %   The input parameters are:
            %       PLAMBDA: reconstructed wavelength
            %       METHOD: mode of doing compressed senssing
            %       KAPPA_VECT_S: the kappa vector matching the current
            %       vector. Usefull to normalized the reconstructed
            %       Polariation

            %
            %   The output parameteres are:
            %       PKAPPA_RECONSTRUCT: the polarization reconstructed and
            %       normalized


              if strcmp(Method,'ContinuousWavelet')==0
                    Pkappa_reconstr = ifftshift(fft(ifftshift(Plambda(1:end-1))));                    
              else
                    %Preal_reconstr = real(ifftshift(fft(ifftshift(Plambda(1:end-1)))));
                    %Pimag_reconstr = imag(ifftshift(fft(ifftshift(Plambda(1:end-1)))));            
                    Preal_reconstr = real(fftshift(ifft(fftshift(Plambda(1:end-1)))));
                    Pimag_reconstr = imag(fftshift(ifft(fftshift(Plambda(1:end-1)))));
              
             end

            % range of kappa in (m^-1)
            conv_factor_kappa = max(kappa_vect_s)-min(kappa_vect_s);

            %Pkappa_reconstr = complex(Preal_reconstr,Pimag_reconstr).*1/conv_factor_kappa;
            Pkappa_reconstr = Pkappa_reconstr.*1/conv_factor_kappa;

         end        
         
          % From wavelength to energy transfer:
         function [DeltaE_grid_unif_J,DeltaE_grid_nonunif_J,lambda_interp_m] = lambda2energ(lambda_vect_m,lambda0_m,K0_m,time_vect_ps)
        %LAMBDA2ENERG(THETA_I_GRID,LAMBDA2_CUT_M,FIG_PLOT)  This function
        % calculates an unifornely spaced energy grid and a non uniformely spaced
        % wavelength grid
        %
        %   The arguments are the following:
        %       LAMBDA_VECT_M: the uniformely spaced grid of positive final
        %       wavelengths
        %       LAMBDA0_M: the initial wavelength in [m]
        %       K0_M: initial wavevector in [1/m]
        %       TIME_VECT_PS: the time vector calculted directly form the current vector. It is used to evaluate the range of energies that we want to achieve.  
        %
        %   The outputs are the following ones:
        %
        %       E_grid_unif_J: uniform gridding in energy transfer [J]
        %       E_grid_nonunif_J: non uniform gridding in energy transfer[J] calculated from lambda_vect_m
        %       lambda_interp_m: non-uniform gridding in wavelength [m]  
        


                load_chess_parameters;


                %% Energy range:
                Delta_time_ps = max(time_vect_ps)-min(time_vect_ps);


                % calculating the energy resolution corresponding to the
                % spin-echo time calculated from the current basis:
                delta_E_J = SE_hbar*1e12/Delta_time_ps;

              

                % calculating the range of energies corresponding to the range of wavelengths:
                E_min_J = SE_h^2/(2*SE_3hemass)*(1/max(lambda_vect_m)^2-1/(lambda0_m)^2);
                E_max_J = SE_h^2/(2*SE_3hemass)*(1/min(lambda_vect_m)^2-1/(lambda0_m)^2);
                
                % calculate the non uniform energy transfer vector from the
                % lambda vector:
                DeltaE_grid_nonunif_J = SE_h^2/(2*SE_3hemass)*(1./lambda_vect_m.^2-1/(lambda0_m)^2);

                res_new = floor((E_max_J-E_min_J)/delta_E_J);

                if res_new> 2e4
                    res_new = 2e4;
                end    

                %% Energy uniform grid:

                % the evenly spaced energy transfer vector:
                DeltaE_grid_unif_J = linspace(E_min_J,E_max_J,res_new);%/C_J2meV;
                E_grid_pos_index=find(DeltaE_grid_unif_J>0);
                Offset=DeltaE_grid_unif_J(E_grid_pos_index(1));
                DeltaE_grid_unif_J=DeltaE_grid_unif_J(2:length(DeltaE_grid_unif_J))-Offset;

                lambda_interp_m = 2*pi./sqrt(2*SE_3hemass*DeltaE_grid_unif_J./SE_hbar^2+K0_m^2);


        end
         
         
          % From wavelength distribution towards scattering function       
         function [SKw_tot,SKw_tot_interp,Jacob_det] = Plambda2Sthetaw_v2(P_lambda,lambda_interp_m,DeltaE_grid_unif_J,E_grid_nonunif_J,lambda0_m,alpha,Method,Coeff,AdditionalInfo)
                %PLAMBDA2SKW(P_LAMBDA,LAMBDA2_POS_M,LAMBDA_2_POS_INTERP_M,THETA_I_GRID,THETA_I_SCAN_M) 
                % In this function, we calculate the scattering function from the
                % wavelength matrix in a uniformely spaced gridding of (DeltaE) for each DeltaK. We multiply the  
                % wavelength distribution by the Jacobian.
                %
                %
                %   The arguments of function are:
                %       P_LAMBDA: the wavelength distribution matrix in terms
                %       (theta_i_interp,lambda2_interp) in [1/m^2*rad^2]
                %       LAMBDA2_POS_M: the positive final wavelengths obtained from the current vector, in [m] (uniformely spaced)
                %       LAMBDA_INTERP_M: the positive interpolated and non-uniformely spaced final wavelengths
                %       [m]
                %       DELTAE_GRID_UNIF_J: the uniformely spaced Delta E gridding in [J]
                %       LAMBDA0_M: incoming wavelength [m]
                %       ALPHA: direction of the line of measurement : I2/I1 = tan(alpha1) in [rad]
                %       COEFF: if the CS reconstructed method is Continuous, this
                %          contains the Fourier coefficients (the amplitude of each
                %            component of the Fourier series).
                %       ADDITIONALINFO: additional info for CS
                %       METHOD: mode of doing compressed senssing
                %
                %   The outputs of the function are:
                %       SKW_TOT: the matrix containing the scattering function in [1/m.J.rad^2]
                %       P_LAMBDA_INTEP: the interpolated wavelength
                %       distribtuion (if the reconstruction method is
                %       'Continuous' the wavelength distribution is
                %       evaluated on a non-uniform wavelength grid, which
                %       corresponds to a unifrom energy grid).  [1/m^2*rad^2]
                %       JACOB_DET: the jacobian  in [m/J]
                %
                    load_chess_parameters;
                    
                    % Calculate the SKW on the non-uniform grid:
                          % The Jacobian to go from Plambda to SKw, in [m/J] 
                       Jacob_det = ...
                            SE_h*lambda0_m^3*SE_3hemass./(2*SE_3hemass*E_grid_nonunif_J*(lambda0_m)^2+(SE_h^2)).^(1.5);


                       % S(theta_i,DeltaE). Note that this is not properly the
                       % scattering function because it depends on theta_i and not on
                       % the momentum transfer! The units are [1/m.J.rad^2]
                       SKw_tot = P_lambda.*Jacob_det';

                
                       % do interpolatin of SKw on a uniform energy grid
                       % (usefull later to calculate the ISF)
                    if strcmp(Method,'ContinuousWavelet')==0
                            SKw_tot_interp = ...
                                    interp1(E_grid_nonunif_J,abs(P_lambda),DeltaE_grid_unif_J,'linear'); 
                    else
                        alpha1 = alpha-pi/2;
                        lambda_interp_back2prefinal=cos(alpha1)*(lambda_interp_m-lambda0_m*tan(alpha1));
                        P_lambda_interp_lambda=Compressed_sensing.EvaluateFunctionFromWaveletCoefficients(lambda_interp_back2prefinal,Coeff,AdditionalInfo);
                        P_lambda_interp_lambda=abs(P_lambda_interp_lambda.');
                        
                        Jacob_det = ...
                            SE_h*lambda0_m^3*SE_3hemass./(2*SE_3hemass*DeltaE_grid_unif_J*(lambda0_m)^2+(SE_h^2)).^(1.5);


                       % S(theta_i,DeltaE). Note that this is not properly the
                       % scattering function because it depends on theta_i and not on
                       % the momentum transfer! The units are [1/m.J.rad^2]
                       SKw_tot_interp = P_lambda_interp_lambda.*Jacob_det;
                        
                    end
                    
                  


         end
         
         % From  scattering function towards intermediate scattering function
         function [ IKt_singleK,time_vect_ps] = Ithetat_calc(Sthetaw_tot,DeltaE_grid_unif_J)
                %IKT_CALC: calculates the intermediate scattering function, i.e. performs the FFT of the S(DeltaK,DeltaE)
                %   The arguments are the following:
                %       SKW_TOT: the scattering function, a matrix in terms of
                %       (deltaK,deltaE) in [1/J.m.rad]
                %       DELTAE_GRID_UNIF_J: the uniformely spaced grid of energy in [J]
                %
                %  The outputs are:     
                %       IKT_TOT: the intermediate scattering function (a matrix in terms of
                %       deltaK and time)
                %       TIME_VECT_ps: the time vector, in [ps]
                %
                %


                    load_chess_parameters;    

                    res_energ = length(DeltaE_grid_unif_J);

                    % limits of the energy range:
                    Energ_max = max(DeltaE_grid_unif_J);% in [J] (max(DeltaE_grid_unif_meV)/conv.C_J2meV/sc.SE_h)*1e-12;% in [ps^{-1}]
                    Energ_min = min(DeltaE_grid_unif_J);% in [J] (min(DeltaE_grid_unif_meV)/conv.C_J2meV/sc.SE_h)*1e-12;% in [ps^{-1}]

                    delta_energ = (Energ_max- Energ_min)/res_energ; % in [J]
                    delta_time=1/(res_energ*delta_energ); % in [1/J]

                    % limits of the time range for the FFT
                    time_min = -round(res_energ/2)*delta_time; % in [1/J]
                    time_max = time_min + (res_energ-1)*delta_time;% in [1/J]

                    % integers:
                    time_min_int = round(time_min/delta_time);
                    time_max_int = round(time_max/delta_time);





                    [IKt_singleK, time_vect_invJ] = Compressed_sensing.OneDimDFTZeroCenteredWithOptions(Sthetaw_tot, Energ_min,Energ_max ,[ time_min_int,time_max_int]);

                     %IKt_singleK = fftshift(ifft(fftshift(Sthetaw_tot)));

                    time_vect_ps = (time_vect_invJ*SE_hbar)*1e12; % in [ps]


         end

         %% Plotting functions:
           % plot the figure for the cs results:
        function [fig_s] = Create_display_figure(fh)
            %CREATE_DISPLAY_FIGURE(figure_h) This function generates the figure on
            %which we will plot the measured and the calculated quantities. This
            %function also sets the pop-up menus for the linear/log scales!
            %
            %   The input is:
            %       FIGURE_H: the figure handle
            %
            %   The output is: 
            %       FIGURE_STRUCT: a structure containing the handles of the axis



                % Create figure
                figure1 = figure(fh);

                % Set the dimensions for the figures (so the pop-up menus lie in the
                % correct position):    
                set(gcf,'Position',[98 130 1120 714]);

                % Create axes for polarization
                axhandle1 = axes('Parent',figure1,...
                    'Position',[0.075 0.53 0.34 0.39]);
                box(axhandle1,'on');
                hold(axhandle1,'all');
                % measured data with positive current:
                plothandle1_1_1=plot(axhandle1,NaN,NaN,'or'); % creates a handle for the  data (in this case the real pol.) to represent in the axis whose handle is axhandle1. 
                plothandle1_2_1=plot(axhandle1,NaN,NaN,'xb'); % Note that it is important to create the handles before labelling the axis since one operation cancels the other
                % measured data with negative current;
                plothandle1_1_2=plot(axhandle1,NaN,NaN,'og'); % creates a handle for the  data (in this case the real pol.) to represent in the axis whose handle is axhandle1. 
                plothandle1_2_2=plot(axhandle1,NaN,NaN,'xk'); 
                % averaged data with positive current:
                plothandle1_1_3=plot(axhandle1,NaN,NaN,'-r');   
                plothandle1_2_3=plot(axhandle1,NaN,NaN,'-b');
                % averaged data with negative current:
                plothandle1_1_4=plot(axhandle1,NaN,NaN,'-g');
                plothandle1_2_4=plot(axhandle1,NaN,NaN,'-k');
                 % if test_mode is on, this handle displays the original polarization data:
                plothandle1_1_5=plot(axhandle1,NaN,NaN,'-r');
                plothandle1_2_5=plot(axhandle1,NaN,NaN,'-b');


                %hold(axhandle1,'all');
                xlabel('Tilted current / A');
                ylabel('Real and Imag Polarisation');
                title('Reconstructed Polarization')


                % Create axes for the wavelength distribution
                axhandle2 = axes('Parent',figure1,...
                    'Position',[0.57 0.54 0.38 0.38]);
                box(axhandle2,'on');
                hold(axhandle2,'on');
                plothandle2_1=plot(axhandle2,NaN,NaN,'ro');
                plothandle2_2=plot(axhandle2,NaN,NaN,'bo');
                xlabel('\lambda_F / m');
                ylabel('\rho(\lambda_F) / m^{-1}');
                title('Final wavelength distribution');


                 % Create axes for the scattering function 
                axhandle3 = axes('Parent',figure1,...
                    'Position',[0.57 0.1 0.38 0.35]);
                box(axhandle3,'on');
                hold(axhandle3,'all');
                plothandle3_1=plot(axhandle3,NaN,NaN,'or'); % real
                plothandle3_2=plot(axhandle3,NaN,NaN,'ob'); % imaginary
                % the wavelength interpolated distribution
                plothandle3_1_2=plot(axhandle3,NaN,NaN,'-k'); % real
                plothandle3_2_2=plot(axhandle3,NaN,NaN,'-k'); % imaginary
               % xlim([-10 30]);
                xlabel('\Delta E [meV]');
                ylabel('S(\Delta E)_\theta in [1/meV*rad]');
                title('Scattering function at a constant \theta');


                % Create axes for the intermediate scattering function
                axhandle4 = axes('Parent',figure1,...
                    'Position',[0.075 0.1 0.35 0.32]);
                hold(axhandle4,'all');
                box(axhandle4,'on');
                plothandle4_1_1=plot(axhandle4,NaN,NaN,'or'); % calculated real IKt
                plothandle4_1_2=plot(axhandle4,NaN,NaN,'xb'); % calculated imaginary IKt
                plothandle4_2_1=plot(axhandle4,NaN,NaN,'-k'); % polarization real
                plothandle4_2_2=plot(axhandle4,NaN,NaN,'-r'); % polarization imaginary
                xlabel('t [ps]');
                ylabel('I(t)_\theta');
                title('Intermediate Scattering function at a constant \theta');



                % store the handlers in a structure

                field1 = 'axhandle1'; value1 = axhandle1;
                field2 = 'axhandle2'; value2 = axhandle2;
                field3 = 'axhandle3'; value3 = axhandle3;
                field4 = 'axhandle4'; value4 = axhandle4;
                field5 = 'plothandle1_1_1'; value5 = plothandle1_1_1;
                field6 = 'plothandle1_2_1'; value6 = plothandle1_2_1;
                field7 = 'plothandle1_1_2'; value7 = plothandle1_1_2;
                field8 = 'plothandle1_2_2'; value8 = plothandle1_2_2;
                field9 = 'plothandle1_1_3'; value9 = plothandle1_1_3;
                field10 = 'plothandle1_2_3'; value10 = plothandle1_2_3;
                field11 = 'plothandle1_1_4'; value11 = plothandle1_1_4;
                field12 = 'plothandle1_2_4'; value12 = plothandle1_2_4;
                field13 = 'plothandle2_1'; value13 = plothandle2_1;
                field14 = 'plothandle2_2'; value14 = plothandle2_2;
                field15 = 'plothandle3_1'; value15 = plothandle3_1;
                field16 = 'plothandle3_2'; value16 = plothandle3_2;
                field17 = 'plothandle3_1_2'; value17 = plothandle3_1_2;
                field18 = 'plothandle3_2_2'; value18 = plothandle3_2_2;
                field19 = 'plothandle4_1_1'; value19 = plothandle4_1_1;
                field20 = 'plothandle4_1_2'; value20 = plothandle4_1_2;
                field21 = 'plothandle4_2_1'; value21 = plothandle4_2_1;
                field22 = 'plothandle4_2_2'; value22 = plothandle4_2_2;
                field23 = 'plothandle1_1_5'; value23 = plothandle1_1_5;
                field24 = 'plothandle1_2_5'; value24 = plothandle1_2_5;
                field25 = 'fh'; value25 = fh;   



                fig_s = struct(field1,value1,field2,value2,field3,value3,field4,value4,field5,value5,field6,value6,...
                    field7,value7,field8,value8,field9,value9,field10,value10,field11,value11,field12,value12,field13,value13,...
                    field14,value14,field15,value15,field16,value16,field17,value17,field18,value18,field19,value19,field20,...
                    value20,field21,value21,field22,value22,field23,value23,field24,value24,field25,value25);

                  %% ===================== Set buttons of the GUI related to log/linear axis: ==============================================


                %  Construct the components for the log/linear popup menu for
                %  polarization
                hpopup1 = uicontrol('Parent', fig_s.fh,'Style','popupmenu',...
                    'String',{'x/y','Linear/Linear','Log/Linear','Linear/Log','Log/Log'},'Position', [75 325 120 20],'Callback',{@popup2_menu_Callback});

                %  Construct the components for the log/linear popup menu for
                % IKT
                hpopup2 = uicontrol('Parent', fig_s.fh,'Style','popupmenu',...
                    'String',{'x/y','Linear/Linear','Log/Linear','Linear/Log','Log/Log'},'Position', [650 30 120 20],'Callback',{@popup3_menu_Callback});

                %  Construct the components for the log/linear popup menu for
                % SKw
                hpopup3 = uicontrol('Parent',fig_s.fh,'Style','popupmenu',...
                    'String',{'x/y','Linear/Linear','Log/Linear','Linear/Log','Log/Log'},'Position', [70 30 120 20],'Callback',{@popup4_menu_Callback});

                %Construct the components for the log/linear popup menu for
                % the wavelength distribution
                hpopup4 = uicontrol('Parent', fig_s.fh,'Style','popupmenu',...
                    'String',{'x/y','Linear/Linear','Log/Linear','Linear/Log','Log/Log'},'Position', [950 325 120 20],'Callback',{@popup5_menu_Callback});


                % show figure toolbar:
                set(fig_s.fh,'toolbar','figure')

                 % Functions which activate the pop-up menus relate to linear and log scale in the axis of the polarization, SF and ISF

                
                function popup1_menu_Callback(source,~,~)
                    % Determine the selected data set.
                    str = get(source, 'String');
                    val = get(source,'Value');

                    % Set current data to the selected data set.
                    switch str{val};
                        case 'Linear/Linear' % User selected linear plot
                            set(fig_s.axhandle1,'XScale','Linear');

                        case 'Log/Linear' % User selects linear/log plot
                            set(fig_s.axhandle1,'XScale','Log');

                        case 'Linear/Log' % User selects linear/log plot
                            set(fig_s.axhandle1,'YScale','Log');  

                        case 'Log/Log' % User selects linear/log plot
                            set(fig_s.axhandle1,'YScale','Log','XScale','Log');      
                    end  

                 end

                function popup2_menu_Callback(source,~,~)

                    % Determine the selected data set.
                    str = get(source, 'String');
                    val = get(source,'Value');

                    % Set current data to the selected data set.
                    switch str{val};
                        case 'Linear/Linear' % User selected linear plot
                            set(fig_s.axhandle1,'XScale','Linear');

                        case 'Log/Linear' % User selects linear/log plot
                            set(fig_s.axhandle1,'XScale','Log');

                        case 'Linear/Log' % User selects linear/log plot
                            set(fig_s.axhandle1,'YScale','Log');  

                        case 'Log/Log' % User selects linear/log plot
                            set(fig_s.axhandle1,'YScale','Log','XScale','Log'); 
                    end

                end

                function popup3_menu_Callback(source,~,~)

                    % Determine the selected data set.
                    str = get(source, 'String');
                    val = get(source,'Value');

                    % Set current data to the selected data set.
                    switch str{val};
                            case 'Linear/Linear' % User selected linear plot
                            set(fig_s.axhandle1,'XScale','Linear');

                        case 'Log/Linear' % User selects linear/log plot
                            set(fig_s.axhandle1,'XScale','Log');

                        case 'Linear/Log' % User selects linear/log plot
                            set(fig_s.axhandle1,'YScale','Log');  

                        case 'Log/Log' % User selects linear/log plot
                            set(fig_s.axhandle1,'YScale','Log','XScale','Log');  
                    end

                end

                function popup4_menu_Callback(source,~,~)

                    % Determine the selected data set.
                    str = get(source, 'String');
                    val = get(source,'Value');

                    % Set current data to the selected data set.
                    switch str{val};
                            case 'Linear/Linear' % User selected linear plot
                            set(fig_s.axhandle1,'XScale','Linear');

                        case 'Log/Linear' % User selects linear/log plot
                            set(fig_s.axhandle1,'XScale','Log');

                        case 'Linear/Log' % User selects linear/log plot
                            set(fig_s.axhandle1,'YScale','Log');  

                        case 'Log/Log' % User selects linear/log plot
                            set(fig_s.axhandle1,'YScale','Log','XScale','Log');  
                    end

                end

        end

        % updates the figure of the cs results:
        function update_display_figure(varargin)
            %UPDATE_DISPLAY_FIGURE(figure_h) This function updates the handles in
            %figure_struct for the next scan.
            %
            %   The input is:
            %      FIGURE_STRUCT: the figure handle
            %      CS_RESULT_STRUCT: structure containing the reconstructed
            %      data to plot. It is provided by the function Pol2ISF
            %      (go there for the detail of the entries)
            %       MEAS: the meas structure of the preexisting data for
            %       comparison with the reconstruction.
            %
            for i=1:2:length(varargin)
                    switch varargin{i}
                        case 'fig_s'
                            fig_s = varargin{i+1};
                        case 'cs_result_struct'
                           cs_result_struct = varargin{i+1};     
                        case 'meas'
                            meas = varargin{i+1};                                                                                      
                        otherwise
                            warning(['variable ' varargin{i} ' is not defined']);
                    end
            end

                if exist('fig_s','var')==0
                    error('No figure structure!');
                end
                
                 if exist('cs_result_struct','var')==0
                    error('No results to plot!');
                 end
                
                 if exist('meas','var')==1
                     test_mode = 1;
                    disp('Test mode is on and the reconstructed Polarization will be compared with the original one');
                end


               %% handles for polarization:

               % Remember that to avoid artifacts in the FFT we remove the
               % last point of Pol_reconstruct in the function Plambda2Pkappa_single_theta
                set(fig_s.plothandle1_1_1,'XData',cs_result_struct.ibasetoreconstruct(1:end-1),'YData',real(cs_result_struct.Pol_reconstr)); % real
                set(fig_s.plothandle1_2_1,'XData',cs_result_struct.ibasetoreconstruct(1:end-1),'YData',imag(cs_result_struct.Pol_reconstr)); % imaginary

                if (test_mode == 1)                    
                    set(fig_s.plothandle1_1_5,'XData',meas.ibase,'YData',meas.mean.Preal); % real
                    set(fig_s.plothandle1_2_5,'XData',meas.ibase,'YData',meas.mean.Pimag); % imaginary
                end
                
                %% handles for wavelength distribution


               set(fig_s.plothandle2_1,'XData',cs_result_struct.lambda_m'*1e10,'YData',real(cs_result_struct.Plambda)'*1e-10); % real
               set(fig_s.plothandle2_2,'XData',cs_result_struct.lambda_m'*1e10,'YData',imag(cs_result_struct.Plambda)'*1e-10); % imaginary

               %set(fig_s.plothandle2_1,'XData',meas.loop(loop_index).lambdaF_m_interp'*1e10,'YData',real(meas.loop(loop_index).Plambda_interp)'*1e-10);
               %set(fig_s.plothandle2_2,'XData',meas.loop(loop_index).lambdaF_m_interp'*1e10,'YData',imag(meas.loop(loop_index).Plambda_interp)'*1e-10);


                %% handles for the scattering function

               set(fig_s.plothandle3_1,'XData',cs_result_struct.Energ_meV,'YData',real(cs_result_struct.SKw)./cs_result_struct.conv.C_J2meV); % real
               set(fig_s.plothandle3_2,'XData',cs_result_struct.Energ_meV,'YData',imag(cs_result_struct.SKw)./cs_result_struct.conv.C_J2meV); % imaginary


               %% handles for the intermediate scattering function

               set(fig_s.plothandle4_1_1,'XData',cs_result_struct.time_ps,'YData',real(cs_result_struct.IKt)); % real 
               set(fig_s.plothandle4_1_2,'XData',cs_result_struct.time_ps,'YData',imag(cs_result_struct.IKt)); % imaginary


          




        end
                      
       
        %% Miscellaneous functions:
        
        function [Pol_s,kappa_vect_s,ibase_s,list_s] = mirror(varargin)
            %MIRROR(ANALYSIS,MEAS) This function checks wether the range of currents is
            %symmetrical and if it isn't, it symmetrizes it..
            %
            %   Input arguments:
            %       ANALYSIS: structure containing the calculated results of the
            %       measurement
            %       MEAS: structure containing the results of the measurements (raw
            %       data)
            %
            %   Output:
            %       POL_S: complex polarization symmetrized
            %       IBASE_S: current vector symmetrized
            %       LIST_S: symmetrized list of currents
             for i=1:2:length(varargin)
                    switch varargin{i}
                        case 'Pol'
                            Pol = varargin{i+1};
                        case 'kappa_vect'
                           kappa_vect = varargin{i+1};     
                        case 'ibase'
                            ibase = varargin{i+1};            
                        case 'list'
                            list = varargin{i+1};                                                     
                        otherwise
                            warning(['variable ' varargin{i} ' is not defined']);
                    end
                end
           


                if min(kappa_vect) >= 0
                    if exist('ibase','var')==1
                        if size(ibase,1)~=1
                            ibase=ibase.';           
                        end               
                        ibase_s = [-fliplr(ibase),ibase(2:length(ibase))];
                    else
                        ibase_s = zeros(1,100);
                    end

                    if size(kappa_vect,1)~=1
                       kappa_vect=kappa_vect.';           
                    end 
                    kappa_vect_s = [-fliplr(kappa_vect),kappa_vect(2:length(kappa_vect))];       

                    if size(Pol,1)~=1
                       Pol=Pol.';           
                    end 
                    Pol_s = complex([fliplr(real(Pol)),real(Pol(2:length(Pol)))],[-fliplr(imag(Pol)),imag(Pol(2:length(Pol)))]); 

                    if size(list,1)~=1
                      list=list.';           
                    end
                    list_s = [fliplr(list),list(2:length(list))];  

                else
                    if exist('ibase','var')==1
                        ibase_s = ibase;
                    else
                        ibase_s = zeros(1,100);
                    end
                    kappa_vect_s = kappa_vect;
                    Pol_s = Pol; 
                    list_s = list;
                end

                % ibase_s = ibase;
               %  kappa_vect_s = kappa_vect;
               %  Pol_s = Pol; 
                % list_s = list;

        end
            
       
        % to calculate the initial conditions of the helium beam:
        function [lambda0_m,K0_m,Energ_i_J,conv_struct] = calc_E0(Energ_0)
                %CALC_E0= This function calculates the initial energy of the beam.
                %
                %   The argument is
                %       ENERG_0: initial energy read on the meas.beam.E0 entry in [meV] 
                % 
                %   The outputs are
                %       LAMBDA0_M: incident wavelength in [m]
                %       K0_M: incident momentum transfer in [1/m]
                %       ENERG_I_J: incident energy in [J]
                %       CONV_STRUCT: structure containing different conversion factors

                    load_chess_parameters;

                    % Conversion constant from Joules to meV:
                    C_J2meV = 6.24e21; % in meV/J


                    % Conversion constant from meV to Hz:
                    C_J2Hz = 1/SE_hbar; % in [Hz.J^-1=1/s*J]

                    % Conversion from J to m^-2:
                    C_J2m2 = SE_hbar^2/(2*SE_3hemass); % in [J.m^2]
                    C_J2rad2m2 = SE_h^2/(2*SE_3hemass); % in [J.rad.m^2]


                    % Conversion from meV to m^1:
                    C_meV2m2 = C_J2m2*C_J2meV; % in [meV.m^2]
                    C_meV2rad2m2 = C_J2rad2m2*C_J2meV; % in [meV.rad.m^2]



                     % Conversion constant from Joules^1 to s:
                    C_invJ2s = SE_hbar; % in J.s


                    % Store the results in a structure
                    conv_struct.C_J2meV = C_J2meV;
                    conv_struct.C_J2Hz = C_J2Hz;
                    conv_struct.C_J2m2 = C_J2m2;
                    conv_struct.C_J2rad2m2 = C_J2rad2m2;
                    conv_struct.C_meV2m2 = C_meV2m2;
                    conv_struct.C_meV2rad2m2 =  C_meV2rad2m2;
                    conv_struct.C_invJ2s =  C_invJ2s;


                    %% Calculate energy, wavelength and wavector:

                    E0 = Energ_0; % Energy of the incoming helium beam in [meV].  

                    K0_Ang = beamprops('energy',E0,3); %wavevector of the incoming beam in [1/m]
                    lambda0_Ang = 2*pi/K0_Ang; %in [Angstroms]
                    lambda0_m = lambda0_Ang*1e-10; % wavelength of an incoming he3atom in [m]. This is in axis lambda1.
                    K0_m = 2*pi/lambda0_m; % the initial momentum in [1/Angstroms]. Warning, it is not projected!!!
                    %V0 = SE_hbar*K0_m/SE_3hemass; % velocity of an incoming he3 atom


                    Energ_i_J = Energ_0/C_J2meV; % in [J]
                   % Energ_i_meV = C_meV2m2*K0_m^2; % in [meV]


        end
                
      

        % average and concatenates measured polarization from different
        % loops:
        function  [meas,list_new] = concat_avg_data(meas,loop_index)
                %CONCAT_DATA this function concatenates the rawdata of the different loops
                %
                %  The input arguments are:
                %       MEAS: the meas structure containing the data for each loop
                %       LOOP_INDEX: index of the current loop
                %
                %
                %  The  ouptut arguments are:
                %       LIST_NEW: updated list of 0s and 1s with the
                %       non-measured/measured currents up to the loop_index
                %       
                %       MEAS: updated version of meas.mean entries
                %           MEAS.MEAN.PREAL
                %           MEAS.MEAN.PIMAG
                %           MEAS.MEAN.Pmag 
                %           MEAS.MEAN.deltaPhase 
                %           MEAS.MEAN.ibase
                %


                    index_list_loop_tot = [];
                    Preal_tot = [];
                    Pimag_tot = [];
                    Pmag_tot = [];
                    deltaPhase_tot = [];
                    ibase_tot = [];

                      %% ===================  Concatenation ===================    
                    for index = 1:loop_index       

                        index_list_loop_tot = [meas.CS_s.list_index_matrix(:,index)' index_list_loop_tot];

                        Preal_tot = [meas.loop(index).Preal,Preal_tot];
                        Pimag_tot = [meas.loop(index).Pimag,Pimag_tot];
                        Pmag_tot = [meas.loop(index).Pmag,Pmag_tot];
                        deltaPhase_tot = [meas.loop(index).deltaPhase,deltaPhase_tot];
                        ibase_tot = [meas.CS_s.ibase_matrix(:,index)' ibase_tot];    
                    end

                   
                    [index_list_loop_sort,iindex] = sort(index_list_loop_tot);

                    Preal_tot = Preal_tot(iindex);
                    Pimag_tot = Pimag_tot(iindex);
                    Pmag_tot = Preal_tot(iindex);
                    deltaPhase_tot = Pimag_tot(iindex);
                    ibase_tot = ibase_tot(iindex);



                    % counting the frequency of each current index (histogram)
                    uniq_index = unique(index_list_loop_sort);
                    count_index = histc(index_list_loop_sort,uniq_index);

                    Preal_avg = zeros(1,length(count_index));
                    Pimag_avg = zeros(1,length(count_index));
                    Pmag_avg = zeros(1,length(count_index));
                    deltaPhase_avg = zeros(1,length(count_index));
                    ibase_avg = zeros(1,length(count_index));




                    % doing the average of the repeated current values 
                    for jj = 1:length(count_index)
                        index_ini = find(index_list_loop_sort == uniq_index(jj),1,'first');
                        Preal_avg(jj) = mean(Preal_tot(index_ini:index_ini+count_index(jj)-1));
                        Pimag_avg(jj) = mean(Pimag_tot(index_ini: index_ini+count_index(jj)-1));
                        Pmag_avg(jj) = mean(Pmag_tot(index_ini:index_ini+count_index(jj)-1));
                        deltaPhase_avg(jj) = mean(deltaPhase_tot(index_ini: index_ini+count_index(jj)-1));
                        ibase_avg(jj) = ibase_tot(index_ini);
                    end

                    [ibase_sort,iibase] = sort(ibase_avg);
                    meas.mean.Preal = Preal_avg(iibase);
                    meas.mean.Pimag = Pimag_avg(iibase);
                    meas.mean.Pmag = Pmag_avg(iibase);
                    meas.mean.deltaPhase = deltaPhase_avg(iibase);
                    meas.mean.ibase = ibase_sort;



                    % prepare the final list for compressed senssing:
                    list_new = zeros(1,length(meas.CS_s.ibasetoreconstruct));
                    list_new(uniq_index) = 1;



        end
         
         % to check the memory storage      
        function [result] = linMem_vI()
            %Copyright (c) 2012, Michael Hirsch
            %All rights reserved.
            %
            %Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
            %
            % Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
            % Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
            %
            %THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

                if ~isunix && ~ismac
                     memStats.error = true;
                    return
                end

                % get the parent process id
                if ismac    
                    [sts,ppid] = unix(['ps -p $PPID -l | ' awkCol('PPID') ]);

                    if sts
                        memStats.error = true;
                        return
                    else %no error

                        % get memory used by the parent process (resident set size)
                        [sts,thisused] = unix(['ps -O rss -p ' strtrim(ppid) ' | awk ''NR>1 {print$2}'' ']);
                        % convert from KB to MB
                        thisused = str2double(thisused)*1024*1e-6; %in MB

                         result = thisused;

                    end


                elseif isunix    
                    [sts,msg] = unix('free -m | grep Mem:');

                    if sts
                        memStats.error = true;
                        return
                    else %no error
                        mems = cell2mat(textscan(msg,'%*s %u %u %u %*u %*u %*u','delimiter',' ','collectoutput',true,'multipleDelimsAsOne',true));
                        memStats.freeMB = mems(3);
                        memStats.usedMB = mems(2);
                        memStats.totalMB = mems(1);
                        memStats.error = false;

                        result = mems(2);
                    end


                end

        end %function
        
        
        
     end
     
end