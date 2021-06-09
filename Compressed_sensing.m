classdef Compressed_sensing
    % FULL_LOOP_FUNC is the library of functions which allow to perform
    % the full loop from the polarization to the intermediate scattering
    % function in the compressed sensing method
    
     properties(Constant)
     end
     
     methods(Static)
         
        % CS with DFT:
         
       function [Recon,ReconDomain_vect]=FastFourierReconFromSamplesZeroCentered(FourierData,SamplingList,Resolution,Fourier_min,Fourier_max)
            % fastfourierreconfromsamples Version 1.0
            % Computes a Compressed Sensing reconstruction from Fourier samples using
            % the FFT.
            % res = resolution of output 
            % kappamin = minimum frequency in terms of kappa
            % kappamax = maximum frequency in terms of kappa
            % lambdamin = Smallest Expected lambda value (if lambdamin and res are not
            %               chosed correctly then there will be wraparound problems)
            % fdata = vector of Fourier samples
            % list = sampling mask 
            row_flag=0;
            if size(FourierData,1)==1
                FourierData=FourierData.';
                row_flag=1;
            end
            if (Resolution~=length(SamplingList))
                error('Size of List Must Match The Resolution');
            end
            if (length(FourierData)~=sum(SamplingList))
                error('List does not correspond to fdata');
            end



            % Compute lambda vector
            Fourier_delta=(Fourier_max-Fourier_min)/(Resolution-1);
            Fourier_max=Fourier_max+Fourier_delta;

            % Compute output domain vector
            ReconDomain_delta=1/(Fourier_delta*Resolution);
            if mod(Resolution,2)==0
                reshalf=Resolution/2;
                ReconDomain_vectpart2=ReconDomain_delta*(0:reshalf).';
                ReconDomain_vectpart1=ReconDomain_delta*(-reshalf+1:-1).';
                ReconDomain_vect=[ReconDomain_vectpart1;ReconDomain_vectpart2];
            else
                reshalf=(Resolution-1)/2;
                ReconDomain_vectpart2=ReconDomain_delta*(0:reshalf).';
                ReconDomain_vectpart1=ReconDomain_delta*(-reshalf:-1).';
                ReconDomain_vect=[ReconDomain_vectpart1;ReconDomain_vectpart2];
            end
            ReconDomain_min=ReconDomain_vect(1);
            ReconDomain_max=ReconDomain_vect(1)+Resolution*ReconDomain_delta;



            % Scale Fourier data before putting into solver 
            fdatascaler=1/(Fourier_max-Fourier_min)*ones(Resolution,1);
            jvect=0:Resolution-1;
            jvect=-2*pi*1i*ReconDomain_min*(Fourier_min+Fourier_delta*jvect);
            jvect=(exp(jvect)).';
            fdatascaler=fdatascaler.*jvect;
            fdatascaler=fdatascaler.^(-1);
            fdatascaler=fdatascaler(SamplingList==1);
            FourierData=FourierData.*fdatascaler;

            % Run Compressed Sensing Solver
            opts = spgSetParms('iterations',100,'decTol',0.01);
            handle=@(input,mode) Compressed_sensing.fastfourier(input,Resolution,SamplingList,mode);
            max_val=max(abs(FourierData));
            FourierData=FourierData/max_val;
            Recon=spg_bpdn(handle,FourierData,0,opts);
            Recon=Recon*max_val;

            % Scale Output Data
            outputscaler=ones(Resolution,1);
            kvect=0:Resolution-1;
            kvect=(-2*pi*1i*Fourier_min*kvect)/(Fourier_max-Fourier_min);
            kvect=(exp(kvect)).';
            outputscaler=outputscaler.*kvect;
            outputscaler=outputscaler.^(-1);
            Recon=Recon.*outputscaler;
            if row_flag==1
                Recon=Recon.';
            end
       end      
       
       function [Recon,ReconDomain_vect]=fastfourierreconfromsamples(FourierData,SamplingList,res,Fourier_min,Fourier_max,ReconDomain_min)
            % fastfourierreconfromsamples Version 1.0
            % Computes a Compressed Sensing reconstruction from Fourier samples using
            % the FFT.
            % res = resolution of output 
            % kappamin = minimum frequency in terms of kappa
            % kappamax = maximum frequency in terms of kappa
            % lambdamin = Smallest Expected lambda value (if lambdamin and res are not
            %               chosed correctly then there will be wraparound problems)
            % fdata = vector of Fourier samples
            % list = sampling mask 
            row_flag=0;
            if size(FourierData,1)==1
                FourierData=FourierData.';
                row_flag=1;
            end
            if (res~=length(SamplingList))
                error('Size of List Must Match The Resolution');
            end
            if (length(FourierData)~=sum(SamplingList))
                error('List does not correspond to fdata');
            end
            % Compute lambda vector
            Fourier_delta=(Fourier_max-Fourier_min)/(res-1);
            Fourier_max=Fourier_max+Fourier_delta;
            ReconDomain_delta=1/(Fourier_delta*res);
            ReconDomain_max=ReconDomain_min+res*ReconDomain_delta;
            ReconDomain_vect=(ReconDomain_min:ReconDomain_delta:ReconDomain_max-ReconDomain_delta)';


            % Scale Fourier data before putting into solver 
            fdatascaler=1/(Fourier_max-Fourier_min)*ones(res,1);
            jvect=0:res-1;
            jvect=-2*pi*1i*ReconDomain_min*(Fourier_min+Fourier_delta*jvect);
            jvect=(exp(jvect)).';
            fdatascaler=fdatascaler.*jvect;
            fdatascaler=fdatascaler.^(-1);
            fdatascaler=fdatascaler(SamplingList==1);
            FourierData=FourierData.*fdatascaler;

            % Run Compressed Sensing Solver
            opts = spgSetParms('iterations',100,'decTol',0.01);
            handle=@(input,mode) Compressed_sensing.fastfourier(input,res,SamplingList,mode);
            max_val=max(abs(FourierData));
            FourierData=FourierData/max_val;
            Recon=spg_bpdn(handle,FourierData,0,opts);
            Recon=Recon*max_val;

            % Scale Output Data
            outputscaler=ones(res,1);
            kvect=0:res-1;
            kvect=(-2*pi*1i*Fourier_min*kvect)/(Fourier_max-Fourier_min);
            kvect=(exp(kvect)).';
            outputscaler=outputscaler.*kvect;
            outputscaler=outputscaler.^(-1);
            Recon=Recon.*outputscaler;
            if row_flag==1
                Recon=Recon.';
            end
        end   

       function output=fastfourier(input,res,list,mode)
        % This is a function that represents the fast fourier transform matrix.
        % input=data to be transformed
        % res=Resolution of final reconstruction
        % list=Samples Taken (see subsampleinlevels.m)
        % mode=Transform mode (normal or conjugate transpose)

            if mode==1
                output=fftn(input);
                output=output(list==1);
            else
                inputpre=zeros(res,1);
                inputpre(list==1)=input;
                input=inputpre;
                output=res*ifftn(input);
            end
       end
         
       function [output,outputdomain_vect]=OneDimDFTZeroCenteredWithOptions(input,inputdomain_min,inputdomain_max,Frequency_Range,Sign)
            % Computes an FT from Fourier samples
            %!!!!!! Output is Zero Centered, not Input!!!!!!!
            % res = resolution of output 
            % Frequency_range = range of frequencies to extract (frequencies indexed by
            % integers)
            % Sign = sign of exponential. leave empty for -2*pi*i*..., 1 for 2*pi*i*...

            
            row_flag=0;
            if size(input,1)==1
                input=input.';
                row_flag=1;
            end
            if nargin<=4
                Sign=-1;
            end

            res=length(input);    
            inputdomain_delta =(inputdomain_max-inputdomain_min)/res;
            % the input domain vector
            inputdomain_vect = inputdomain_min+(0:res-1)*inputdomain_delta;  
            % shift the input for suitably doing the fft 
            %input = [input(inputdomain_vect>=inputdomain_vect(find(input==max(input))));input(inputdomain_vect<inputdomain_vect(find(input==max(input))))];
            % Compute output domain vector
            outputdomain_delta=1/(inputdomain_delta*res);
            if mod(res,2)==0
                reshalf=res/2;
                outputdomain_vectpart2=outputdomain_delta*(0:reshalf).';
                outputdomain_vectpart1=outputdomain_delta*(-reshalf+1:-1).';
                outputdomain_vect=[outputdomain_vectpart1;outputdomain_vectpart2];
            else
                reshalf=(res-1)/2;
                outputdomain_vectpart2=outputdomain_delta*(0:reshalf).';
                outputdomain_vectpart1=outputdomain_delta*(-reshalf:-1).';
                outputdomain_vect=[outputdomain_vectpart1;outputdomain_vectpart2];
            end
            outputdomain_min=outputdomain_vect(1);
            outputdomain_max=outputdomain_vect(1)+res*outputdomain_delta;

            % Scale and Shift Output Data
            inputscaler=ones(res,1);
            kvect=0:res-1;
            kvect=(Sign*2*pi*1i*outputdomain_min*kvect)/(outputdomain_max-outputdomain_min);
            kvect=(exp(kvect)).';
            inputscaler=inputscaler.*kvect;
            input=input.*inputscaler;


            % Run FT
            max_val=max(abs(input));
            input=input/max_val;
            if Sign==1
                output=res*ifftn(input);
            else
                output=fftn(input);
            end
            output=output*max_val;


            % Scale and Shift Fourier data
            outputscaler=1/(outputdomain_max-outputdomain_min)*ones(res,1);
            jvect=0:res-1;
            jvect=Sign*2*pi*1i*inputdomain_min*(outputdomain_min+outputdomain_delta*jvect);
            jvect=(exp(jvect)).';
            outputscaler=outputscaler.*jvect;
            output=output.*outputscaler;

            if nargin>3
                frequency_vect=round(outputdomain_vect/outputdomain_delta);
                frequency_select=(frequency_vect>=Frequency_Range(1))&(frequency_vect<=Frequency_Range(2));
                outputdomain_vect=outputdomain_vect(frequency_select);
                output=output(frequency_select);
            end

            if row_flag==1
                output=output.';
            end
       end 
       
       function Output=BlockFFT(Input,FFTRes,OutputRes,Sign)
            % Computes FFT for an inputs/outputs that don't fit a specific DFT size.
            % This allows the DFT size to be specified separately.

            % Fill Input with zeros so that the size is a multiple of FFTRes
            InputFill=FFTRes-mod(length(Input),FFTRes);
            Input=[Input;zeros(InputFill,1)];

            InputMatrix=reshape(Input,FFTRes,[]);

            % fft applies 1D DFT to each column separately here
            if FFTRes~=1
                if Sign==-1
                    FFTMatrix=fft(InputMatrix);
                else
                    FFTMatrix=ifft(InputMatrix)*FFTRes;
                end
                FFTVect=sum(FFTMatrix,2);
            else
                FFTVect=sum(InputMatrix,2);
            end

            % Repeat Output and select entries as necessary so that it has the same
            % size as OutputRes
            NumberOfCopies=ceil(OutputRes/FFTRes);
            FFTVect=repmat(FFTVect,NumberOfCopies,1);
            Output=FFTVect(1:OutputRes);
       end 
        
       function wcoef = FWT_PO(x,L,qmf)
            % FWT_PO -- Forward Wavelet Transform (periodized, orthogonal)
            %  Usage
            %    wc = FWT_PO(x,L,qmf)
            %  Inputs
            %    x    1-d signal; length(x) = 2^J
            %    L    Coarsest Level of V_0;  L << J
            %    qmf  quadrature mirror filter (orthonormal)
            %  Outputs
            %    wc    1-d wavelet transform of x.
            %
            %  Description
            %    1. qmf filter may be obtained from MakeONFilter   
            %    2. usually, length(qmf) < 2^(L+1)
            %    3. To reconstruct use IWT_PO
            %
            %  See Also
            %    IWT_PO, MakeONFilter
            %
            n = length(x) ;
            J=1;
            while abs(n*2^(-J)-ceil(n*2^(-J)))<1e-5
               J=J+1; 
            end
            J=J-1;
            if (J<L)
                error('Cannot divide by 2 enough times');
            end

            aqmf=-( (-1).^(1:length(qmf)) ).*qmf;

            wcoef = zeros(1,n) ;
            beta = x(:).';  %take samples at finest scale as beta-coeffts
            for j=0:1:L-1
                alfa = Compressed_sensing.iconv( aqmf,[ beta( 2:length(beta) ) beta(1) ]);
                m = length(alfa);
                alfa = alfa(1:2:(m-1));

                wcoef((n*2^(-j-1)+1):(n*2^(-j))) = alfa;

                beta = Compressed_sensing.aconv(qmf,beta);
                m = length(beta);
                beta = beta(1:2:(m-1));
            end
            wcoef(1:(n*2^(-L))) = beta;
            wcoef = reshape(wcoef,[size(x,1),size(x,2)]);
        end
       
        % CS with DWT (wavelets basis):
       function output=fastfourierwavelet(input,res,filter,scalenum,list,mode)
          % This is a function that represents the fast fourier transform matrix
        % combined with a discrete wavelet transform matrix.
        % input=data to be transformed
        % res=Resolution of final reconstruction
        % list=Samples Taken (see subsampleinlevels.m)
        % mode=Transform mode (normal or conjugate transpose)
        % filter=wavelet filter
        % scalenum=number of wavelet levels
            inputshift=0;
            shiftfound=0;
            dyadscale=2^scalenum;
            while shiftfound==0
               if abs(floor((res-inputshift)/dyadscale)-(res-inputshift)/dyadscale)<1e-10
                   shiftfound=1;
               else
                   inputshift=inputshift+1;
               end
            end
            if mode==1
                inputpre=input(1:inputshift);
                inputchop=input(inputshift+1:length(input));
                inputchop=IWT_PO(inputchop,scalenum,filter);
                input=[inputpre;inputchop];
                output=fftn(input);
                output=output(list==1);
            else
                inputpre=zeros(res,1);
                inputpre(list==1)=input;
                input=inputpre;
                input=res*ifftn(input);
                inputpre=input(1:inputshift);
                inputchop=input(inputshift+1:length(input));
                inputchop=Compressed_sensing.FWT_PO(inputchop,scalenum,filter);
                output=[inputpre;inputchop];
            end

        end
            
       function Output=FastFourierWithWavelets(Input,Resolution,SamplingList,Filter,Levels,Mode)
            % This is a function that represents the fast fourier transform matrix.
            % input=data to be transformed
            % res=Resolution of final reconstruction
            % list=Samples Taken (see subsampleinlevels.m)
            % mode=Transform mode (normal or conjugate transpose)

            InputShift=0;
            ShiftFound=0;
            Dyadscale=2^Levels;
            while ShiftFound==0
               if abs(floor((Resolution-InputShift)/Dyadscale)-(Resolution-InputShift)/Dyadscale)<1e-10
                   ShiftFound=1;
               else
                   InputShift=InputShift+1;
               end
            end
            if Mode==1
                InputPre=Input(1:InputShift);
                InputChop=Input(InputShift+1:length(Input));
                InputChop=Compressed_sensing.IWT_PO(InputChop,Levels,Filter);
                Input=[InputPre;InputChop];
                Output=fftn(Input);
                Output=Output(SamplingList==1);
            else
                InputPre=zeros(Resolution,1);
                InputPre(SamplingList==1)=Input;
                Input=InputPre;
                Input=Resolution*ifftn(Input);
                InputPre=Input(1:InputShift);
                InputChop=Input(InputShift+1:length(Input));
                InputChop=Compressed_sensing.FWT_PO(InputChop,Levels,Filter);
                Output=[InputPre;InputChop];
            end
        end
        
       function [Recon,ReconDomain_vect]=FastFourierReconFromSamplesZeroCenteredWithWavelets(FourierData,SamplingList,Resolution,Fourier_min,Fourier_max,Filter,Levels)
            % fastfourierreconfromsamples Version 1.0
            % Computes a Compressed Sensing reconstruction from Fourier samples using
            % the FFT.
            % res = resolution of output 
            % kappamin = minimum frequency in terms of kappa
            % kappamax = maximum frequency in terms of kappa
            % lambdamin = Smallest Expected lambda value (if lambdamin and res are not
            %               chosed correctly then there will be wraparound problems)
            % fdata = vector of Fourier samples
            % list = sampling mask 
            
            row_flag=0;
            if size(FourierData,1)==1
                FourierData=FourierData.';
                row_flag=1;
            end
            if (Resolution~=length(SamplingList))
                error('Size of List Must Match The Resolution');
            end
            if (length(FourierData)~=sum(SamplingList))
                error('List does not correspond to fdata');
            end



            % Compute lambda vector
            Fourier_delta=(Fourier_max-Fourier_min)/(Resolution-1);
            Fourier_max=Fourier_max+Fourier_delta;

            % Compute output domain vector
            ReconDomain_delta=1/(Fourier_delta*Resolution);
            if mod(Resolution,2)==0
                reshalf=Resolution/2;
                ReconDomain_vectpart2=ReconDomain_delta*(0:reshalf).';
                ReconDomain_vectpart1=ReconDomain_delta*(-reshalf+1:-1).';
                ReconDomain_vect=[ReconDomain_vectpart1;ReconDomain_vectpart2];
            else
                reshalf=(Resolution-1)/2;
                ReconDomain_vectpart2=ReconDomain_delta*(0:reshalf).';
                ReconDomain_vectpart1=ReconDomain_delta*(-reshalf:-1).';
                ReconDomain_vect=[ReconDomain_vectpart1;ReconDomain_vectpart2];
            end
            ReconDomain_min=ReconDomain_vect(1);
            ReconDomain_max=ReconDomain_vect(1)+Resolution*ReconDomain_delta;



            % Scale Fourier data before putting into solver 
            fdatascaler=1/(Fourier_max-Fourier_min)*ones(Resolution,1);
            jvect=0:Resolution-1;
            jvect=-2*pi*1i*ReconDomain_min*(Fourier_min+Fourier_delta*jvect);
            jvect=(exp(jvect)).';
            fdatascaler=fdatascaler.*jvect;
            fdatascaler=fdatascaler.^(-1);
            fdatascaler=fdatascaler(SamplingList==1);
            FourierData=FourierData.*fdatascaler;

            % Run Compressed Sensing Solver
            opts = spgSetParms('iterations',100,'decTol',0.01);
            handle=@(Input,Mode) Compressed_sensing.FastFourierWithWavelets(Input,Resolution,SamplingList,Filter,Levels,Mode);
            max_val=max(abs(FourierData));
            FourierData=FourierData/max_val;
            Recon=spg_bpdn(handle,FourierData,0,opts);
            Recon=Recon*max_val;
            ReconShift=0;
            ShiftFound=0;
            Dyadscale=2^Levels;
            while ShiftFound==0
               if abs(floor((Resolution-ReconShift)/Dyadscale)-(Resolution-ReconShift)/Dyadscale)<1e-10
                   ShiftFound=1;
               else
                   ReconShift=ReconShift+1;
               end
            end
            ReconPre=Recon(1:ReconShift);
            ReconChop=Recon(ReconShift+1:length(Recon));
            ReconChop=IWT_PO(ReconChop,Levels,Filter);
            Recon=[ReconPre;ReconChop];

            % Scale Output Data
            outputscaler=ones(Resolution,1);
            kvect=0:Resolution-1;
            kvect=(-2*pi*1i*Fourier_min*kvect)/(Fourier_max-Fourier_min);
            kvect=(exp(kvect)).';
            outputscaler=outputscaler.*kvect;
            outputscaler=outputscaler.^(-1);
            Recon=Recon.*outputscaler;
            if row_flag==1
                Recon=Recon.';
            end
        end
        
       function Output=WaveletFourierDiagonal(FrequencyMin,FrequencyMax,Filter,Iterations,LowestScale,HighestScale)
    % function in the compressed sensing method
                %  Computes the fourier transform of the scaling function and mother wavelet
                %  evaluated at 2^(-j)m where m=FrequencyMin,...,FrequencyMax
                %  and j=0,...,HighestScale.
                %  The Fourier Transform uses the -2*pi*1i convention.
                %  The Scaling function and wavelet are supported on [0,Length(Filter)-1]
                %  Iterations refers to number of iterative steps
                
                        Support=length(Filter); Degree=Support/2;

                        % Changed to rephrase the filter in so that it now refers to the scaling
                        % function phi(x) expressed in terms of its scaled version 2*phi(2x) and its
                        % shifts (i.e. not  2^{1/2}*phi(2x))
                        Filter=Filter/sqrt(2);

                        FrequencyTotal=FrequencyMax-FrequencyMin+1;
                        Output=zeros(FrequencyTotal,HighestScale-LowestScale+1);
                        SupportRun=(0:Support-1);
                        FrequencyVect=(FrequencyMin:FrequencyMax).';
                        Exponent=-2*pi*1i*FrequencyVect;

                        % Scaling function term
                        NewVect=ones(FrequencyTotal,1);
                        for k=1:Iterations
                            IterExponent=2^(-k-LowestScale)*Exponent*SupportRun;
                            IterMatrix=exp(IterExponent);
                            IterVect=IterMatrix*Filter.';
                            NewVect=NewVect.*IterVect;
                        end
                        Output(:,1)=NewVect;

                        % Wavelet Terms
                        for j=LowestScale:HighestScale
                            ExponentScale=0.5*Exponent.*2^(-j);
                            ExponentFilter=2*pi*1i*(0.5*2^(-j)*FrequencyVect + 0.5);

                            % The Degree term here is to shift the wavelet back to have support on
                            % [0,Length(Filter)-1] since the ExponentialFilter terms flips the
                            % support of the scaling function around 0.
                            MatrixFilter=ExponentFilter*(SupportRun-2*(Degree-1));
                            MatrixFilter=exp(MatrixFilter);
                            NewVect=MatrixFilter*Filter';
                            for k=1:Iterations
                                IterMatrix=2^(-k)*ExponentScale*SupportRun;
                                IterMatrix=exp(IterMatrix);
                                IterMatrix=IterMatrix*Filter';
                                NewVect=NewVect.*IterMatrix;
                            end
                            NewVect=exp(ExponentScale).*NewVect;
                            Output(:,2+j-LowestScale)=NewVect;
                        end


        end
                
       function Output=PeriodicWavelet2FourierHandle(Input,LowestScale,HighestScale,FrequencyMin,FrequencyMax,List,BasisFData,Mode)
                    FrequencyTotal=FrequencyMax-FrequencyMin+1;
                    if Mode==1 %Wavelet Coefficients to Fourier Coefficients
                        Input=Compressed_sensing.WaveletVector2LevelsForm(Input,LowestScale,HighestScale);
                        Output=zeros(FrequencyTotal,1);

                        % Base Scaling Function Component
                        ScalingFData=BasisFData(:,1);
                        j=LowestScale;
                        % Where the FFT work is done
                        % Strictly Negative Frequencies
                        DFTOutput1=Compressed_sensing.BlockFFT(Input{1},2^j,-FrequencyMin+1,1);
                        DFTOutput1=DFTOutput1(2:length(DFTOutput1));
                        % Positive Frequencies
                        DFTOutput2=Compressed_sensing.BlockFFT(Input{1},2^j,FrequencyMax+1,-1);
                        DFTOutput=[flipud(DFTOutput1);DFTOutput2];

                        Output=Output+2^(-j/2)*ScalingFData.*DFTOutput;

                        % Wavelet Levels
                        for j=LowestScale:HighestScale
                            WaveletFData=BasisFData(:,j-LowestScale+2);

                            % Where the FFT work is done
                            % Strictly Negative Frequencies
                            DFTOutput1=Compressed_sensing.BlockFFT(Input{j-LowestScale+2},2^j,-FrequencyMin+1,1);
                            DFTOutput1=DFTOutput1(2:length(DFTOutput1));
                            % Positive Frequencies

                            DFTOutput2=Compressed_sensing.BlockFFT(Input{j-LowestScale+2},2^j,FrequencyMax+1,-1);
                            DFTOutput=[flipud(DFTOutput1);DFTOutput2];

                            Output=Output+2^(-j/2)*WaveletFData.*DFTOutput;
                        end
                        Output=Output(List==1);
                    else %Fourier Coefficients to Wavelet Coefficients
                        InputFull=zeros(FrequencyTotal,1);
                        InputFull(List==1)=Input;
                        Input=InputFull;

                        Output=cell(HighestScale-LowestScale+2,1);

                        j=LowestScale;
                        % Base Scaling Function Component
                        ScalingFData=conj(BasisFData(:,1));
                        OutputPreFT=Input.*ScalingFData;

                        OutputPreFT1=OutputPreFT(1:-FrequencyMin+1);
                        OutputPreFT1(-FrequencyMin+1)=0;
                        DFTOutput1=Compressed_sensing.BlockFFT(flipud(OutputPreFT1),2^j,2^j,-1);

                        OutputPreFT2=OutputPreFT(-FrequencyMin+1:FrequencyTotal);
                        DFTOutput2=Compressed_sensing.BlockFFT(OutputPreFT2,2^j,2^j,1);

                        % Where the FFT work is done 
                        DFTOutput=DFTOutput1+DFTOutput2;

                        Output{1}=2^(-j/2)*DFTOutput;

                        % Wavelet Levels
                        for j=LowestScale:HighestScale
                            WaveletFData=conj(BasisFData(:,j-LowestScale+2));
                            OutputPreFT=Input.*WaveletFData;

                            OutputPreFT1=OutputPreFT(1:-FrequencyMin+1);
                            OutputPreFT1(-FrequencyMin+1)=0;
                            DFTOutput1=Compressed_sensing.BlockFFT(flipud(OutputPreFT1),2^j,2^j,-1);

                            OutputPreFT2=OutputPreFT(-FrequencyMin+1:FrequencyTotal);
                            DFTOutput2=Compressed_sensing.BlockFFT(OutputPreFT2,2^j,2^j,1);

                            % Where the FFT work is done 
                            DFTOutput=DFTOutput1+DFTOutput2;

                            Output{j-LowestScale+2}=2^(-j/2)*DFTOutput;
                        end
                        Output=Compressed_sensing.WaveletLevels2VectorForm(Output,LowestScale,HighestScale);
                    end
        end

       function Output=WaveletLevels2VectorForm(Input,LowestScale,HighestScale)
            Output=zeros(2^(HighestScale+1),1);
            Output(1:2^LowestScale)=Input{1};
            for j=LowestScale:HighestScale
               Output(2^j+1:2^(j+1))=Input{j-LowestScale+2};
            end 
       end

       function Output=WaveletVector2LevelsForm(Input,LowestScale,HighestScale)
            Output=cell(HighestScale-LowestScale+2,1);
            Output{1}=Input(1:2^LowestScale);
            for j=LowestScale:HighestScale
                Output{j-LowestScale+2}=Input(2^j+1:2^(j+1));
            end
       end
        
       function x = IWT_PO(wc,L,qmf)
        % IWT_PO -- Inverse Wavelet Transform (periodized, orthogonal)
        %  Usage
        %    x = IWT_PO(wc,L,qmf)
        %  Inputs
        %    wc     1-d wavelet transform: length(wc) = 2^J.
        %    L      Coarsest scale (2^(-L) = scale of V_0); L << J;
        %    qmf    quadrature mirror filter
        %  Outputs
        %    x      1-d signal reconstructed from wc
        %
        %  Description
        %    Suppose wc = FWT_PO(x,L,qmf) where qmf is an orthonormal quad. mirror
        %    filter, e.g. one made by MakeONFilter. Then x can be reconstructed by
        %      x = IWT_PO(wc,L,qmf)
        %
        %  See Also
        %    FWT_PO, MakeONFilter
        %
                n = length(wc) ;
                J=1;
                while abs(n*2^(-J)-ceil(n*2^(-J)))<1e-5
                   J=J+1; 
                end
                J=J-1;
                if (J<L)
                    error('Cannot divide by 2 enough times');
                end

                aqmf=-( (-1).^(1:length(qmf)) ).*qmf;

                wcoef  = wc(:).';
                x= wcoef(1:n*2^(-L));
                for j=0:L-1
                    x1=Compressed_sensing.iconv( qmf, Compressed_sensing.UpSampleN(x) );
                    xj=wcoef((n*2^(j-L)+1):(n*2^(j-L+1)));
                    xj= Compressed_sensing.UpSampleN(xj);
                    m = length(xj);
                    xj = [ xj(m) xj( 1: (m-1) )];
                    x2=Compressed_sensing.aconv(aqmf, xj );
                    x=x1+x2;
                end

                x = reshape(x,[size(wc,1),size(wc,2)]);

                %
                % Copyright (c) 1993. Iain M. Johnstone
                %     
        end
        
        % Continuous CS: 
        function Output=EvaluateFunctionFromWaveletCoefficients(Input,Coeff,AdditionalInfo,ImagingScale)
                Filter=AdditionalInfo.Filter;
                DomainRange=AdditionalInfo.DomainRange;
                LowestScale=AdditionalInfo.LowestScale;
                HighestScale=AdditionalInfo.HighestScale;
                if nargin<=3
                   ImagingScale=HighestScale+5; 
                end
                small=(DomainRange(2)-DomainRange(1))/1000;
                DomainCheck=(Input>DomainRange(2)+small) | (Input<DomainRange(1)-small);
                % display(AdditionalInfo.DomainRange)
                % display([min(Input),max(Input)])
                if sum(DomainCheck)>=1
                    error('Requested Arguments For the Function Fall out of Reconstruction Domain');
                end

                Input=(Input-DomainRange(1))/(DomainRange(2)-DomainRange(1));

                FilterShift=length(Filter)/2-1;
                Output=zeros(2^ImagingScale,1);
                Output(1:2^(HighestScale+1))=Coeff;

                OutputLevel=Output(1:2^LowestScale);
                Output(1:2^LowestScale)=OutputLevel;
                for j=LowestScale:HighestScale
                    OutputLevel=Output(2^j+1:2^(j+1));
                    OutputLevel=circshift(OutputLevel,FilterShift);
                    Output(2^j+1:2^(j+1))=OutputLevel;
                end
                Output=Compressed_sensing.IWT_PO(Output,ImagingScale-LowestScale,Filter);
                Output=2^(ImagingScale/2)*Output;

                Input=max(round(Input*2^ImagingScale),1);
                Output=Output(Input);

        end

        function [Coefficients,AdditionalInfo]=...
            SolveForWaveletCoefficientsFromFourierSamples(FourierData,SamplingList,Filter,FrequencyRange,FrequencyStep,ReconDomainMin,LowestScale,HighestScale,Sign)
             % Solves for Wavelet Coefficients From Fourier Samples    
            FrequencyTotal=FrequencyRange(2)-FrequencyRange(1)+1;
            if nargin<=7
                HighestScale=ceil((log(FrequencyTotal)/log(2)))+2;
            end
            if nargin<=6
                LowestScale=0;
            end
            if nargin<=8
                Sign=1;
            end

            row_flag=0;
            if size(FourierData,1)==1
                FourierData=FourierData.';
                row_flag=1;
            end
            if (length(FourierData)~=sum(SamplingList))
                error('List does not correspond to fdata');
            end


            if Sign==-1
                FrequencyRange=fliplr(-FrequencyRange);
                SamplingList=fliplr(SamplingList);
                FourierData=fliplr(FourierData);
            end

            FrequencyMin=FrequencyRange(1); FrequencyMax=FrequencyRange(2);

            Frequency_Vect=((FrequencyMin:FrequencyMax)*FrequencyStep).';
            FrequencyShift=2*pi*1i*Frequency_Vect(SamplingList==1)*ReconDomainMin;
            FourierData=FrequencyStep*FourierData.*exp(FrequencyShift);

            % Run Compressed Sensing Solver
            BasisFData=Compressed_sensing.WaveletFourierDiagonal(FrequencyMin,FrequencyMax,Filter,20,LowestScale,HighestScale);
            SolverOptions = spgSetParms('iterations',100,'decTol',0.01);
            MatrixHandle=@(Input,Mode) Compressed_sensing.PeriodicWavelet2FourierHandle(Input,LowestScale,HighestScale,FrequencyMin,FrequencyMax,SamplingList,BasisFData,Mode);

            max_val=max(abs(FourierData));
            FourierData=FourierData/max_val;
            Recon=spg_bpdn(MatrixHandle,FourierData,0,SolverOptions);
            Recon=Recon*max_val;

            if row_flag==1
                Recon=Recon.';
            end
            Coefficients=Recon;
            DomainRange=[0,1/FrequencyStep]+ReconDomainMin;
            AdditionalInfo.DomainRange=DomainRange;
            AdditionalInfo.Filter=Filter;
            AdditionalInfo.LowestScale=LowestScale;
            AdditionalInfo.HighestScale=HighestScale;
        end
        
               
        % Miscellaneous functions:
            

        function y = iconv(f,x)
            % iconv -- Convolution Tool for Two-Scale Transform
            %  Usage
            %    y = iconv(f,x)
            %  Inputs
            %    f   filter
            %    x   1-d signal
            %  Outputs
            %    y   filtered result
            %
            %  Description
            %    Filtering by periodic convolution of x with f

            n = length(x);
            p = length(f);
            if p <= n,
               xpadded = [x((n+1-p):n) x];
            else
               z = zeros(1,p);
               for i=1:p,
                   imod = 1 + rem(p*n -p + i-1,n);
                   z(i) = x(imod);
               end
               xpadded = [z x];
            end
            ypadded = filter(f,1,xpadded);
            y = ypadded((p+1):(n+p));
        end

        function y = aconv(f,x)
                % aconv -- Convolution Tool for Two-Scale Transform
                %  Usage
                %    y = aconv(f,x)
                %  Inputs
                %    f    filter
                %    x    1-d signal
                %  Outputs
                %    y    filtered result
                %
                %  Description
                %    Filtering by periodic convolution of x with the
                %    time-reverse of f.
                %
                %  See Also
                %    iconv, UpDyadHi, UpDyadLo, DownDyadHi, DownDyadLo
                %

                n = length(x);
                p = length(f);
                if p < n,
                   xpadded = [x x(1:p)];
                else
                   z = zeros(1,p);
                   for i=1:p,
                       imod = 1 + rem(i-1,n);
                       z(i) = x(imod);
                   end
                   xpadded = [x z];
                end
                fflip = f(length(f):-1:1);
                ypadded = filter(fflip,1,xpadded);
                y = ypadded(p:(n+p-1));
        end
                
        function y = UpSampleN(x,s)
                % UpSample -- Upsampling operator
                %  Usage
                %    u = UpSample(d[,s]) 
                %  Inputs
                %    d   1-d signal, of length n
                %    s   upsampling scale, default = 2
                %  Outputs
                %    u   1-d signal, of length s*n with zeros
                %        interpolating alternate samples
                %        u(s*i-1) = d(i), i=1,...,n
                %
                if nargin == 1, s = 2; 
                end
                n = length(x)*s;
                y = zeros(1,n);
                y(1:s:(n-s+1) )=x;
        end                  
            
       function f = MakeONFilter(Type,Par)
            % MakeONFilter -- Generate Orthonormal QMF Filter for Wavelet Transform
            %  Usage
            %    qmf = MakeONFilter(Type,Par)
            %  Inputs
            %    Type   string, 'Haar', 'Beylkin', 'Coiflet', 'Daubechies',
            %           'Symmlet', 'Vaidyanathan','Battle'
            %    Par    integer, it is a parameter related to the support and vanishing
            %           moments of the wavelets, explained below for each wavelet.
            %
            % Outputs
            %    qmf    quadrature mirror filter
            %
            %  Description
            %    The Haar filter (which could be considered a Daubechies-2) was the
            %    first wavelet, though not called as such, and is discontinuous.
            %
            %    The Beylkin filter places roots for the frequency response function
            %    close to the Nyquist frequency on the real axis.
            %
            %    The Coiflet filters are designed to give both the mother and father
            %    wavelets 2*Par vanishing moments; here Par may be one of 1,2,3,4 or 5.
            %
            %    The Daubechies filters are minimal phase filters that generate wavelets
            %    which have a minimal support for a given number of vanishing moments.
            %    They are indexed by their length, Par, which may be one of
            %    4,6,8,10,12,14,16,18 or 20. The number of vanishing moments is par/2.
            %
            %    Symmlets are also wavelets within a minimum size support for a given 
            %    number of vanishing moments, but they are as symmetrical as possible,
            %    as opposed to the Daubechies filters which are highly asymmetrical.
            %    They are indexed by Par, which specifies the number of vanishing
            %    moments and is equal to half the size of the support. It ranges 
            %    from 4 to 10.
            %
            %    The Vaidyanathan filter gives an exact reconstruction, but does not
            %    satisfy any moment condition.  The filter has been optimized for
            %    speech coding.
            %
            %    The Battle-Lemarie filter generate spline orthogonal wavelet basis.
            %    The parameter Par gives the degree of the spline. The number of 
            %    vanishing moments is Par+1.
            %
            %  See Also
            %    FWT_PO, IWT_PO, FWT2_PO, IWT2_PO, WPAnalysis
            %
            %  References
            %    The books by Daubechies and Wickerhauser.
            %

            if strcmp(Type,'Haar'),
                f = [1 1] ./ sqrt(2);
            end

            if strcmp(Type,'Beylkin'),
                f = [	.099305765374	.424215360813	.699825214057	...
                        .449718251149	-.110927598348	-.264497231446	...
                        .026900308804	.155538731877	-.017520746267	...
                        -.088543630623	.019679866044	.042916387274	...
                        -.017460408696	-.014365807969	.010040411845	...
                        .001484234782	-.002736031626	.000640485329	];
            end

            if strcmp(Type,'Coiflet'),
                if Par==1,
                    f = [	.038580777748	-.126969125396	-.077161555496	...
                            .607491641386	.745687558934	.226584265197	];
                end
                       % MakeONFilter -- Generate Orthonormal QMF Filter for Wavelet Transform
            %  Usage
            %    qmf = MakeONFilter(Type,Par)
            %  Inputs
            %    Type   string, 'Haar', 'Beylkin', 'Coiflet', 'Daubechies',
            %           'Symmlet', 'Vaidyanathan','Battle'
            %    Par    integer, it is a parameter related to the support and vanishing
            %           moments of the wavelets, explained below for each wavelet.
            %
            % Outputs
            %    qmf    quadrature mirror filter
            %
            %  Description
            %    The Haar filter (which could be considered a Daubechies-2) was the
            %    first wavelet, though not called as such, and is discontinuous.
            %
            %    The Beylkin filter places roots for the frequency response function
            %    close to the Nyquist frequency on the real axis.
            %
            %    The Coiflet filters are designed to give both the mother and father
            %    wavelets 2*Par vanishing moments; here Par may be one of 1,2,3,4 or 5.
            %
            %    The Daubechies filters are minimal phase filters that generate wavelets
            %    which have a minimal support for a given number of vanishing moments.
            %    They are indexed by their length, Par, which may be one of
            %    4,6,8,10,12,14,16,18 or 20. The number of vanishing moments is par/2.
            %
            %    Symmlets are also wavelets within a minimum size support for a given 
            %    number of vanishing moments, but they are as symmetrical as possible,
            %    as opposed to the Daubechies filters which are highly asymmetrical.
            %    They are indexed by Par, which specifies the number of vanishing
            %    moments and is equal to half the size of the support. It ranges 
            %    from 4 to 10.
            %
            %    The Vaidyanathan filter gives an exact reconstruction, but does not
            %    satisfy any moment condition.  The filter has been optimized for
            %    speech coding.
            %
            %    The Battle-Lemarie filter generate spline orthogonal wavelet basis.
            %    The parameter Par gives the degree of the spline. The number of 
            %    vanishing moments is Par+1.
            %
            %  See Also
            %    FWT_PO, IWT_PO, FWT2_PO, IWT2_PO, WPAnalysis
            %
            %  References
            %    The books by Daubechies and Wickerhauser.
            %

            if strcmp(Type,'Haar'),
                f = [1 1] ./ sqrt(2);
            end

            if strcmp(Type,'Beylkin'),
                f = [	.099305765374	.424215360813	.699825214057	...
                        .449718251149	-.110927598348	-.264497231446	...
                        .026900308804	.155538731877	-.017520746267	...
                        -.088543630623	.019679866044	.042916387274	...
                        -.017460408696	-.014365807969	.010040411845	...
                        .001484234782	-.002736031626	.000640485329	];
            end

            if strcmp(Type,'Coiflet'),
                if Par==1,
                    f = [	.038580777748	-.126969125396	-.077161555496	...
                            .607491641386	.745687558934	.226584265197	];
                end
                if Par==2,
                    f = [	.016387336463	-.041464936782	-.067372554722	...
                            .386110066823	.812723635450	.417005184424	...
                            -.076488599078	-.059434418646	.023680171947	...
                            .005611434819	-.001823208871	-.000720549445	];
                end
                if Par==3,
                    f = [	-.003793512864	.007782596426	.023452696142	...
                            -.065771911281	-.061123390003	.405176902410	...
                            .793777222626	.428483476378	-.071799821619	...
                            -.082301927106	.034555027573	.015880544864	...
                            -.009007976137	-.002574517688	.001117518771	...
                            .000466216960	-.000070983303	-.000034599773	];
                end
                if Par==4,
                    f = [	.000892313668	-.001629492013	-.007346166328	...
                            .016068943964	.026682300156	-.081266699680	...
                            -.056077313316	.415308407030	.782238930920	...
                            .434386056491	-.066627474263	-.096220442034	...
                            .039334427123	.025082261845	-.015211731527	...
                            -.005658286686	.003751436157	.001266561929	...
                            -.000589020757	-.000259974552	.000062339034	...
                            .000031229876	-.000003259680	-.000001784985	];
                end
                if Par==5,
                    f = [	-.000212080863	.000358589677	.002178236305	...
                            -.004159358782	-.010131117538	.023408156762	...
                            .028168029062	-.091920010549	-.052043163216	...
                            .421566206729	.774289603740	.437991626228	...
                            -.062035963906	-.105574208706	.041289208741	...
                            .032683574283	-.019761779012	-.009164231153	...
                            .006764185419	.002433373209	-.001662863769	...
                            -.000638131296	.000302259520	.000140541149	...
                            -.000041340484	-.000021315014	.000003734597	...
                            .000002063806	-.000000167408	-.000000095158	];
                end
            end

            if strcmp(Type,'Daubechies'),
                if Par==4,  
                    f = [	.482962913145	.836516303738	...
                            .224143868042	-.129409522551	];
                end
                if Par==6, 
                    f = [	.332670552950	.806891509311	...
                            .459877502118	-.135011020010	...
                            -.085441273882	.035226291882	];
                end
                if Par==8,
                    f = [ 	.230377813309	.714846570553	...
                            .630880767930	-.027983769417	...
                            -.187034811719	.030841381836	...
                            .032883011667	-.010597401785	];
                end
                if Par==10,
                    f = [	.160102397974	.603829269797	.724308528438	...
                            .138428145901	-.242294887066	-.032244869585	...
                            .077571493840	-.006241490213	-.012580751999	...
                            .003335725285									];
                end
                if Par==12,
                    f = [	.111540743350	.494623890398	.751133908021	...
                            .315250351709	-.226264693965	-.129766867567	...
                            .097501605587	.027522865530	-.031582039317	...
                            .000553842201	.004777257511	-.001077301085	];
                end
                if Par==14,
                    f = [	.077852054085	.396539319482	.729132090846	...
                            .469782287405	-.143906003929	-.224036184994	...
                            .071309219267	.080612609151	-.038029936935	...
                            -.016574541631	.012550998556	.000429577973	...
                            -.001801640704	.000353713800					];
                end
                if Par==16,
                    f = [	.054415842243	.312871590914	.675630736297	...
                            .585354683654	-.015829105256	-.284015542962	...
                            .000472484574	.128747426620	-.017369301002	...
                            -.044088253931	.013981027917	.008746094047	...
                            -.004870352993	-.000391740373	.000675449406	...
                            -.000117476784									];
                end
                if Par==18,
                    f = [	.038077947364	.243834674613	.604823123690	...
                            .657288078051	.133197385825	-.293273783279	...
                            -.096840783223	.148540749338	.030725681479	...
                            -.067632829061	.000250947115	.022361662124	...
                            -.004723204758	-.004281503682	.001847646883	...
                            .000230385764	-.000251963189	.000039347320	];
                end
                if Par==20,
                    f = [	.026670057901	.188176800078	.527201188932	...
                            .688459039454	.281172343661	-.249846424327	...
                            -.195946274377	.127369340336	.093057364604	...
                            -.071394147166	-.029457536822	.033212674059	...
                            .003606553567	-.010733175483	.001395351747	...
                            .001992405295	-.000685856695	-.000116466855	...
                            .000093588670	-.000013264203					];
                end
            end

            if strcmp(Type,'Symmlet'),
                if Par==4,
                    f = [	-.107148901418	-.041910965125	.703739068656	...
                            1.136658243408	.421234534204	-.140317624179	...
                            -.017824701442	.045570345896					];
                end
                if Par==5,
                    f = [	.038654795955	.041746864422	-.055344186117	...
                            .281990696854	1.023052966894	.896581648380	...
                            .023478923136	-.247951362613	-.029842499869	...
                            .027632152958									];
                end
                if Par==6,  
                    f = [	.021784700327	.004936612372	-.166863215412	...
                            -.068323121587	.694457972958	1.113892783926	...
                            .477904371333	-.102724969862	-.029783751299	...
                            .063250562660	.002499922093	-.011031867509	];
                end
                if Par==7,
                    f = [	.003792658534	-.001481225915	-.017870431651	...
                            .043155452582	.096014767936	-.070078291222	...
                            .024665659489	.758162601964	1.085782709814	...
                            .408183939725	-.198056706807	-.152463871896	...
                            .005671342686	.014521394762					];
                end
                if Par==8, 
                    f = [	.002672793393	-.000428394300	-.021145686528	...
                            .005386388754	.069490465911	-.038493521263	...
                            -.073462508761	.515398670374	1.099106630537	...
                            .680745347190	-.086653615406	-.202648655286	...
                            .010758611751	.044823623042	-.000766690896	... 
                            -.004783458512									];
                end
                if Par==9,
                    f = [	.001512487309	-.000669141509	-.014515578553	...
                            .012528896242	.087791251554	-.025786445930	...
                            -.270893783503	.049882830959	.873048407349	...
                            1.015259790832	.337658923602	-.077172161097	...
                            .000825140929	.042744433602	-.016303351226	...
                            -.018769396836	.000876502539	.001981193736	];
                end
                if Par==10,
                    f = [	.001089170447	.000135245020	-.012220642630	...
                            -.002072363923	.064950924579	.016418869426	...
                            -.225558972234	-.100240215031	.667071338154	...
                            1.088251530500	.542813011213	-.050256540092	...
                            -.045240772218	.070703567550	.008152816799	...
                            -.028786231926	-.001137535314	.006495728375	...
                            .000080661204	-.000649589896					];
                end
            end

            if strcmp(Type,'Vaidyanathan'),
                f = [	-.000062906118	.000343631905	-.000453956620	...
                        -.000944897136	.002843834547	.000708137504	...
                        -.008839103409	.003153847056	.019687215010	...
                        -.014853448005	-.035470398607	.038742619293	...
                        .055892523691	-.077709750902	-.083928884366	...
                        .131971661417	.135084227129	-.194450471766	...
                        -.263494802488	.201612161775	.635601059872	...
                        .572797793211	.250184129505	.045799334111		];
            end

            if strcmp(Type,'Battle'),
                if Par == 1,
                       g = [0.578163    0.280931   -0.0488618   -0.0367309 ...
                            0.012003    0.00706442 -0.00274588 -0.00155701 ...
                            0.000652922 0.000361781 -0.000158601 -0.0000867523
                    ];
                end

                if Par == 3,

                g = [0.541736    0.30683    -0.035498    -0.0778079 ...
                         0.0226846   0.0297468     -0.0121455 -0.0127154 ...
                         0.00614143 0.00579932    -0.00307863 -0.00274529 ...
                         0.00154624 0.00133086 -0.000780468 -0.00065562 ...
                     0.000395946 0.000326749 -0.000201818 -0.000164264 ...
                         0.000103307
                    ];
                end

                if Par == 5,
                g = [0.528374    0.312869    -0.0261771   -0.0914068 ...
                         0.0208414    0.0433544 -0.0148537 -0.0229951  ...
                         0.00990635 0.0128754    -0.00639886 -0.00746848 ...
                         0.00407882 0.00444002 -0.00258816    -0.00268646 ...
                         0.00164132 0.00164659 -0.00104207 -0.00101912 ...
                    0.000662836 0.000635563 -0.000422485 -0.000398759 ...
                    0.000269842 0.000251419 -0.000172685 -0.000159168 ...
                    0.000110709 0.000101113
                    ];
                end
                    l = length(g);
                    f = zeros(1,2*l-1);
                    f(l:2*l-1) = g;
                    f(1:l-1) = reverse(g(2:l));
            end

            f = f ./ norm(f);

            % 
            % Copyright (c) 1993-5. Jonathan Buckheit and David Donoho
            %     




            %
            %  Part of Wavelab Version 850
            %  Built Tue Jan  3 13:20:40 EST 2006
            %  This is Copyrighted Material
            %  For Copying permissions see COPYING.m
            %  Comments? e-mail wavelab@stat.stanford.edu 
     if Par==2,
                    f = [	.016387336463	-.041464936782	-.067372554722	...
                            .386110066823	.812723635450	.417005184424	...
                            -.076488599078	-.059434418646	.023680171947	...
                            .005611434819	-.001823208871	-.000720549445	];
                end
                if Par==3,
                    f = [	-.003793512864	.007782596426	.023452696142	...
                            -.065771911281	-.061123390003	.405176902410	...
                            .793777222626	.428483476378	-.071799821619	...
                            -.082301927106	.034555027573	.015880544864	...
                            -.009007976137	-.002574517688	.001117518771	...
                            .000466216960	-.000070983303	-.000034599773	];
                end
                if Par==4,
                    f = [	.000892313668	-.001629492013	-.007346166328	...
                            .016068943964	.026682300156	-.081266699680	...
                            -.056077313316	.415308407030	.782238930920	...
                            .434386056491	-.066627474263	-.096220442034	...
                            .039334427123	.025082261845	-.015211731527	...
                            -.005658286686	.003751436157	.001266561929	...
                            -.000589020757	-.000259974552	.000062339034	...
                            .000031229876	-.000003259680	-.000001784985	];
                end
                if Par==5,
                    f = [	-.000212080863	.000358589677	.002178236305	...
                            -.004159358782	-.010131117538	.023408156762	...
                            .028168029062	-.091920010549	-.052043163216	...
                            .421566206729	.774289603740	.437991626228	...
                            -.062035963906	-.105574208706	.041289208741	...
                            .032683574283	-.019761779012	-.009164231153	...
                            .006764185419	.002433373209	-.001662863769	...
                            -.000638131296	.000302259520	.000140541149	...
                            -.000041340484	-.000021315014	.000003734597	...
                            .000002063806	-.000000167408	-.000000095158	];
                end
            end

            if strcmp(Type,'Daubechies'),
                if Par==4,  
                    f = [	.482962913145	.836516303738	...
                            .224143868042	-.129409522551	];
                end
                if Par==6, 
                    f = [	.332670552950	.806891509311	...
                            .459877502118	-.135011020010	...
                            -.085441273882	.035226291882	];
                end
                if Par==8,
                    f = [ 	.230377813309	.714846570553	...
                            .630880767930	-.027983769417	...
                            -.187034811719	.030841381836	...
                            .032883011667	-.010597401785	];
                end
                if Par==10,
                    f = [	.160102397974	.603829269797	.724308528438	...
                            .138428145901	-.242294887066	-.032244869585	...
                            .077571493840	-.006241490213	-.012580751999	...
                            .003335725285									];
                end
                if Par==12,
                    f = [	.111540743350	.494623890398	.751133908021	...
                            .315250351709	-.226264693965	-.129766867567	...
                            .097501605587	.027522865530	-.031582039317	...
                            .000553842201	.004777257511	-.001077301085	];
                end
                if Par==14,
                    f = [	.077852054085	.396539319482	.729132090846	...
                            .469782287405	-.143906003929	-.224036184994	...
                            .071309219267	.080612609151	-.038029936935	...
                            -.016574541631	.012550998556	.000429577973	...
                            -.001801640704	.000353713800					];
                end
                if Par==16,
                    f = [	.054415842243	.312871590914	.675630736297	...
                            .585354683654	-.015829105256	-.284015542962	...
                            .000472484574	.128747426620	-.017369301002	...
                            -.044088253931	.013981027917	.008746094047	...
                            -.004870352993	-.000391740373	.000675449406	...
                            -.000117476784									];
                end
                if Par==18,
                    f = [	.038077947364	.243834674613	.604823123690	...
                            .657288078051	.133197385825	-.293273783279	...
                            -.096840783223	.148540749338	.030725681479	...
                            -.067632829061	.000250947115	.022361662124	...
                            -.004723204758	-.004281503682	.001847646883	...
                            .000230385764	-.000251963189	.000039347320	];
                end
                if Par==20,
                    f = [	.026670057901	.188176800078	.527201188932	...
                            .688459039454	.281172343661	-.249846424327	...
                            -.195946274377	.127369340336	.093057364604	...
                            -.071394147166	-.029457536822	.033212674059	...
                            .003606553567	-.010733175483	.001395351747	...
                            .001992405295	-.000685856695	-.000116466855	...
                            .000093588670	-.000013264203					];
                end
            end

            if strcmp(Type,'Symmlet'),
                if Par==4,
                    f = [	-.107148901418	-.041910965125	.703739068656	...
                            1.136658243408	.421234534204	-.140317624179	...
                            -.017824701442	.045570345896					];
                end
                if Par==5,
                    f = [	.038654795955	.041746864422	-.055344186117	...
                            .281990696854	1.023052966894	.896581648380	...
                            .023478923136	-.247951362613	-.029842499869	...
                            .027632152958									];
                end
                if Par==6,  
                    f = [	.021784700327	.004936612372	-.166863215412	...
                            -.068323121587	.694457972958	1.113892783926	...
                            .477904371333	-.102724969862	-.029783751299	...
                            .063250562660	.002499922093	-.011031867509	];
                end
                if Par==7,
                    f = [	.003792658534	-.001481225915	-.017870431651	...
                            .043155452582	.096014767936	-.070078291222	...
                            .024665659489	.758162601964	1.085782709814	...
                            .408183939725	-.198056706807	-.152463871896	...
                            .005671342686	.014521394762					];
                end
                if Par==8, 
                    f = [	.002672793393	-.000428394300	-.021145686528	...
                            .005386388754	.069490465911	-.038493521263	...
                            -.073462508761	.515398670374	1.099106630537	...
                            .680745347190	-.086653615406	-.202648655286	...
                            .010758611751	.044823623042	-.000766690896	... 
                            -.004783458512									];
                end
                if Par==9,
                    f = [	.001512487309	-.000669141509	-.014515578553	...
                            .012528896242	.087791251554	-.025786445930	...
                            -.270893783503	.049882830959	.873048407349	...
                            1.015259790832	.337658923602	-.077172161097	...
                            .000825140929	.042744433602	-.016303351226	...
                            -.018769396836	.000876502539	.001981193736	];
                end
                if Par==10,
                    f = [	.001089170447	.000135245020	-.012220642630	...
                            -.002072363923	.064950924579	.016418869426	...
                            -.225558972234	-.100240215031	.667071338154	...
                            1.088251530500	.542813011213	-.050256540092	...
                            -.045240772218	.070703567550	.008152816799	...
                            -.028786231926	-.001137535314	.006495728375	...
                            .000080661204	-.000649589896					];
                end
            end

            if strcmp(Type,'Vaidyanathan'),
                f = [	-.000062906118	.000343631905	-.000453956620	...
                        -.000944897136	.002843834547	.000708137504	...
                        -.008839103409	.003153847056	.019687215010	...
                        -.014853448005	-.035470398607	.038742619293	...
                        .055892523691	-.077709750902	-.083928884366	...
                        .131971661417	.135084227129	-.194450471766	...
                        -.263494802488	.201612161775	.635601059872	...
                        .572797793211	.250184129505	.045799334111		];
            end

            if strcmp(Type,'Battle'),
                if Par == 1,
                       g = [0.578163    0.280931   -0.0488618   -0.0367309 ...
                            0.012003    0.00706442 -0.00274588 -0.00155701 ...
                            0.000652922 0.000361781 -0.000158601 -0.0000867523
                    ];
                end

                if Par == 3,

                g = [0.541736    0.30683    -0.035498    -0.0778079 ...
                         0.0226846   0.0297468     -0.0121455 -0.0127154 ...
                         0.00614143 0.00579932    -0.00307863 -0.00274529 ...
                         0.00154624 0.00133086 -0.000780468 -0.00065562 ...
                     0.000395946 0.000326749 -0.000201818 -0.000164264 ...
                         0.000103307
                    ];
                end

                if Par == 5,
                g = [0.528374    0.312869    -0.0261771   -0.0914068 ...
                         0.0208414    0.0433544 -0.0148537 -0.0229951  ...
                         0.00990635 0.0128754    -0.00639886 -0.00746848 ...
                         0.00407882 0.00444002 -0.00258816    -0.00268646 ...
                         0.00164132 0.00164659 -0.00104207 -0.00101912 ...
                    0.000662836 0.000635563 -0.000422485 -0.000398759 ...
                    0.000269842 0.000251419 -0.000172685 -0.000159168 ...
                    0.000110709 0.000101113
                    ];
                end
                    l = length(g);
                    f = zeros(1,2*l-1);
                    f(l:2*l-1) = g;
                    f(1:l-1) = reverse(g(2:l));
            end

            f = f ./ norm(f);

            % 
            % Copyright (c) 1993-5. Jonathan Buckheit and David Donoho
            %     




            %
            %  Part of Wavelab Version 850
            %  Built Tue Jan  3 13:20:40 EST 2006
            %  This is Copyrighted Material
            %  For Copying permissions see COPYING.m
            %  Comments? e-mail wavelab@stat.stanford.edu 

        end
             
     end
end