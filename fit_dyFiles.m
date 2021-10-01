classdef fit_dyFiles
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(Constant)
        projFold='Cu111_4';
    end
    
    methods(Static)
        
        function plot_funcional(x,ffun,gof,modes)
            % The assumpsion made in this function is that the 'modes' are
            % named such that a center gaussian(lorentzian) is named
            % cg_<whatever> (cl_<whatever>). Otherwise, special cases are
            % dealt with.
            cn = coeffnames(ffun);
            cv = coeffvalues(ffun);
            modes = modes(~cellfun('isempty',modes));
            
            hold on; plot(x,ffun(x),'k')
            
            i=1;
            for j=1:length(modes)
                current_mode = modes{j};
                switch current_mode(1:3)
                    case 'cl_'
                        y = (cv(i)/pi)*((cv(i+1)/2)./(x.^2+(cv(i+1)/2)^2));
                        modes{j}= [modes{j} ': x_0=0' sprintf(', sig=%0.4f, A=%0.4f',cv(i+1),cv(i))];
                        i = i+2;
                    case 'ncl'
                        y = (cv(i)/pi)*((cv(i+2)/2)./((x-cv(i+1)).^2+(cv(i+2)/2)^2));
                        modes{j}= [modes{j} ': ' sprintf('x_0=%0.2f, sig=%0.2f, A=%0.4f',cv(i+1),cv(i+2),cv(i))];
                        i = i+3;
                    case 'cg_'
                        y = cv(i)*exp(-1*x.^2/cv(i+1)^2);
                        modes{j}= [modes{j} ': x_0=0, ' sprintf('sig=%0.2f',cv(i+1))];
                        i=i+2;
                    case 'ncg'
                        y = cv(i)*exp(-1*(x-cv(i+1)).^2/cv(i+2)^2);
                        modes{j}= [modes{j} ': ' sprintf('x_0=%0.2f, sig=%0.2f, A=%0.4f',cv(i+1),cv(i+2),cv(i))];
                        i=i+3;
                    case 'lin'
                        modes{j}=[]; continue;
                    case 'sp1' % special
                        for k=0:11, eval([cn{i+k} '=' num2str(cv(i+k)) ';']); end
                        %for k=0:11, eval(['n' num2str(1+k) '=cv(' num2str(i+k) ');']); end
                        y1 = (n1/pi)*((n2/2)./((x-1*n3).^2+(n2/2)^2));
                        y2 = (n1/pi)*((n2/2)./((x+1*n3).^2+(n2/2)^2));
                        y3 = (n4*n1/pi)*((n5*n2/2)./((x-(2*n3+n10)).^2+(n5*n2/2)^2));
                        y4 = (n4*n1/pi)*((n5*n2/2)./((x+(2*n3+n10)).^2+(n5*n2/2)^2));
                        y5 = (n6*n4*n1/pi)*((n7*n5*n2/2)./((x-(3*n3+n11)).^2+(n7*n5*n2/2)^2));
                        y6 = (n6*n4*n1/pi)*((n7*n5*n2/2)./((x+(3*n3+n11)).^2+(n7*n5*n2/2)^2));
                        y7 = (n8*n6*n4*n1/pi)*((n9*n7*n5*n2/2)./((x-(4*n3+n12)).^2+(n9*n7*n5*n2/2)^2));
                        y8 = (n8*n6*n4*n1/pi)*((n9*n7*n5*n2/2)./((x+(4*n3+n12)).^2+(n9*n7*n5*n2/2)^2));
                        y = y1+y2+y3+y4+y5+y6+y7+y8;
                        
                        if n1>=0.01, modes{j}= [modes{j} ': ' sprintf('x_0=%0.2f,width=%0.2f',n3,n2)]; end
                        i=i+12;
                        if cv(i-12)<0.01, modes{j}=[]; continue; end
                        
                    case 'sp2' % special
                        for k=0:8, eval(['n' num2str(1+k) '=cv(' num2str(i+k) ');']); end
                        y1 = (n1/pi)*((n2/2)./((x-1*n3).^2+(n2/2)^2));
                        y2 = (n1/pi)*((n2/2)./((x+1*n3).^2+(n2/2)^2));
                        y3 = (n4*n1/pi)*((n2/2)./((x-(2*n3+n7)).^2+(n2/2)^2));
                        y4 = (n4*n1/pi)*((n2/2)./((x+(2*n3+n7)).^2+(n2/2)^2));
                        y5 = (n5*n4*n1/pi)*((n2/2)./((x-(3*n3+n8)).^2+(n2/2)^2));
                        y6 = (n5*n4*n1/pi)*((n2/2)./((x+(3*n3+n8)).^2+(n2/2)^2));
                        y7 = (n6*n5*n4*n1/pi)*((n2/2)./((x-(4*n3+n9)).^2+(n2/2)^2));
                        y8 = (n6*n5*n4*n1/pi)*((n2/2)./((x+(4*n3+n9)).^2+(n2/2)^2));
                        y = y1+y2+y3+y4+y5+y6+y7+y8;
                        
                        if n1>=0.01, modes{j}= [modes{j} ': ' sprintf('x_0=%0.2f,width=%0.2f',n3,n2)]; end
                        i=i+9;
                        if cv(i-9)<0.01, modes{j}=[]; continue; end
                        
                    case 'sp3' % special
                        for k=0:4, eval(['n' num2str(1+k) '=cv(' num2str(i+k) ');']); end
                        y1 = (n1/pi)*((n2/2)./((x-1*n3).^2+(n2/2)^2));
                        y2 = (n1/pi)*((n2/2)./((x+1*n3).^2+(n2/2)^2));
                        y3 = (n4/pi)*((n5/2)./((x-2*n3).^2+(n5/2)^2));
                        y4 = (n4/pi)*((n5/2)./((x+2*n3).^2+(n5/2)^2));
                        y = y1+y2+y3+y4;
                        
                        modes{j}= [modes{j} ': ' sprintf('x_0=%0.2f,width=%0.2f',n3,n2)];
                        i=i+5;
                        
                    case 'pl4' % a fourth order polynomial
                        for k=0:4, eval(['n' num2str(1+k) '=cv(' num2str(i+k) ');']); end
                        y=n5+n1*x+n2*x.^2+n3*x.^3+n4*x.^4;
                        
                        modes{j}= [modes{j} ': ' sprintf('n_5=%0.2f,n_4=%0.2f',n5,n4)];
                        i=i+5;
                    otherwise
                        y = cv(i)*exp(-1*(x-cv(i+1)).^2/cv(i+2)^2)+cv(i+3)*exp(-1*(x-cv(i+4)).^2/cv(i+5)^2);
                        if cv(i)>0.01, modes{j}= [modes{j} ': ' sprintf('%0.2f',cv(i+1))]; end
                        if cv(i+3)>0.005, modes{j}= [modes{j} ': ' sprintf('%0.2f',cv(i+4))]; end
                        i=i+6;
                        if cv(i-6)<0.01 && cv(i-3)<0.01, modes{j}=[]; continue; end
                end
                hold on; plot(x,y);
            end
            
            modes = modes(~cellfun('isempty',modes));
            legend('data', ['total fit, Rsquare=' sprintf('%0.3f',gof.rsquare) ', rmse=' sprintf('%0.3f',gof.rmse) ], modes{:})
            
        end
        
        function res = prepare_measured_set_for_fitting(varargin)
            
            if nargout < 1, error('no res output param was defined, EXITING'); end
            
            [filenames, data_path] = postprocess_dyfiles.get_files(varargin);
            load_chess_parameters
            
            %% Load measurements
            resInd = 1;
            for i=1:length(filenames)
                
                filename=[data_path char(filenames(i)) '.mat'];
                
                variableInfo = who('-file', filename);
                if ismember('processed_meas', variableInfo)
                    load(filename,'processed_meas')
                    meas = processed_meas;
                else
                    disp(['File ' filename ' lacks the variable processed_meas'])
                    continue;
                end
                
                meas.filename=filename;
                
                % Calculate the standard deviation of each point in the ISF
                % (for the Pmag, and Preal)
                tmp = reshape([meas.loop.Preal],meas.numloops,[]);
                meas.mean.realErr=std(tmp,0,1);
                tmp = reshape([meas.loop.Pmag],meas.numloops,[]);
                meas.mean.magErr=std(tmp,0,1);
                
                res(resInd)=meas;
                resInd = resInd+1;
            end
            
            % Sort by dK, is is more comfortable
            [~,indx]=sort(abs([res.dK]));
            res = res(indx);
            
            if isfield(res, 'endStatus')
                for j=1:length(res)
                    gammaOffSpec = res(j).endStatus.gamma - res(j).endStatus.gammaSpecular;
                    res(j).gammaTrue = gammaOffSpec + 22.2;
                end
            end
            
            if isfield(res, 'beam')
                for j=1:length(res)
                    res(j).E0 = res(j).beam.E0;
                    [ki,~,~] = beamprops('energy',res(j).E0,3);
                    [res(j).dK,~]=k_xfer(ki,res(j).gammaTrue,44.4);
                end
            end
            
        end
        
        function res = prepare_md_set_for_fitting_forJohnsCode(filename) % for john's code
            
            if nargin < 1, error('no filename was given'); end
            if nargout < 1, error('no res output param was defined, EXITING'); end
            
            load_chess_parameters
            load(filename)
            
            if ~exist('omegashift','var')
                omegashift(1:ntscat/2)=omega(ntscat/2+1:ntscat)-ntscat*(omega(2)-omega(1));
                omegashift(ntscat/2+1:ntscat)=omega(1:ntscat/2);
            end
            
            % smooth the data
            tsmoothfwhm=500e1;
            tsmoothsigma=tsmoothfwhm/sqrt(8*log(2));
            smoothvec=exp(-min(time.^2,(time(ntscat)-time).^2)/2/tsmoothsigma^2);
            for i=1:size(intensity,2)
                for j=1:size(intensity,3)
                    res(i,j).filename = filename;
                    
                    max_intensity = max(intensity(:,i,j));
                    tmpSKw = fftshift(fft(ifft(squeeze(intensity(:,i,j)))'.*smoothvec));
                    res(i,j).SKw = tmpSKw/max(tmpSKw)*max_intensity;
                    res(i,j).ISF0 = ISF0(:,i,j);
                    res(i,j).ISF0 = res(i,j).ISF0(1:length(ISF0(:,i,j))/2);
                    res(i,j).ISF0 = res(i,j).ISF0./res(i,j).ISF0(1);
                    res(i,j).setime = time(1:length(ISF0(:,i,j))/2);
                    res(i,j).mean.Preal = res(i,j).ISF0;
                    res(i,j).mean.Pimag = zeros(size(res(i,j).ISF0));
                    res(i,j).Energ_meV = omegashift * SE_hbar / SE_e * 1E3 * 1E12 ;
                    res(i,j).gammaTrue = [];
                    res(i,j).E0 = [];
                    res(i,j).KxKy = DeltaK(:,i,j);
                    % Determine azimuth
                    tmp = res(i,j).KxKy(1)/res(i,j).KxKy(2);
                    if tmp == 0 || isinf(tmp), res(i,j).alpha=110;
                    elseif abs(tmp-sqrt(3))<1e-3, res(i,j).alpha=112;
                    elseif isnan(tmp)
                        res(i,j).alpha=NaN;
                    else
                        res(i,j).alpha=[];
                    end
                    res(i,j).dK = sqrt(DeltaK(1,i,j).*DeltaK(1,i,j)+DeltaK(2,i,j).*DeltaK(2,i,j));
                    res(i,j).temperature = temp;
                end
            end
            
            res = reshape(res,1,[]);
            
            % Sort by dK, is is more comfortable
            [tmp_dK,indx]=sort(abs([res.dK]));
            res = res(indx);
            
        end
        
        function res = prepare_md_set_for_fitting_forSimulink(filename, varargin)
            
            if nargout < 1, error('no res output param was defined, EXITING'); end
            
            % parse varargin
            prsdArgs = inputParser;   % Create instance of inputParser class.
            prsdArgs.addParameter('isf_str','isf_inc_CoM', @ischar);
            prsdArgs.addParameter('decim_fctr',1, @isnumeric);
            prsdArgs.parse(varargin{:});
            isf_str = prsdArgs.Results.isf_str;
            decim_fctr = round(prsdArgs.Results.decim_fctr);
            
            load_chess_parameters
            load(filename,'params')
            
            max_isf_str = ['max_' isf_str];
            
            tmp = load(filename,isf_str);
            tmp1=fieldnames(tmp);
            isf = tmp.(tmp1{1});
            
            tmp = load(filename,max_isf_str);
            if isempty(tmp) || isempty(fieldnames(tmp))
                max_isf = ones(size(isf,1),1,size(isf,3));
            else
                tmp1=fieldnames(tmp);
                max_isf = tmp.(tmp1{1});
            end
            
            % Correct for previous versions
            if ~isfield(params,'dK') && isfield(params,'dk')
                params.dK = params.dk;
            end
            
            % Assume ISF was truncated
            tmp = 0;
            if isfield(params,'t_isf')
                tmp = length(params.t_isf)-length(params.t)/2;
            else
                params.t_isf = params.t;
            end
            if tmp > 1, error('this is unexpected'); end
            isf_cyclic = [isf fliplr(isf(:,1+tmp:end-tmp,:))];
            
            ntscat = size(isf_cyclic,2);
            
            t = 0:params.isf_sample_time:params.isf_sample_time*params.N_ISF_steps-params.thermalizing_time; % total sim time (without thermalization)
            
            omega_max  = 2*pi/params.isf_sample_time;
            dOmega     = 2*pi/t(end);
            omega      = 0:dOmega:omega_max-dOmega;
            omegashift = -omega_max/2:dOmega:omega_max/2-dOmega;
            
            %             omega=linspace(0,2*pi/t(end)*(ntscat-1),ntscat);
            %             omegashift(1:ntscat/2)=omega(ntscat/2+1:ntscat)-ntscat*(omega(2)-omega(1));
            %             omegashift(ntscat/2+1:ntscat)=omega(1:ntscat/2);
            
            zeroIndx =floor(length(omegashift)/2)+1;
            indx_tmp = 1:decim_fctr:zeroIndx;
            indx = [1:decim_fctr:length(omegashift)]+(zeroIndx-indx_tmp(end));
            Energ_meV = omegashift(indx) * SE_hbar / SE_e * 1E3 * 1E12 ;
            
            % Assume the isf is not symmetric (only 'positive' times are provided)
            
            intensity = fft(isf_cyclic,[],2);
            
            % smooth the data
            tsmoothfwhm=500e1;
            tsmoothsigma=tsmoothfwhm/sqrt(8*log(2));
            smoothvec=exp(-min(t.^2,(t(ntscat)-t).^2)/2/tsmoothsigma^2);
            smoothvec = ones(size(smoothvec)); %disp('Skw is not smoothed')
            
            for i=1:size(isf,1) % dK's
                for j=1:size(isf,3) % alpha's
                    res(i,j).filename = filename;
                    max_intensity = max(intensity(i,:,j));
                    tmpSKw = fftshift(fft(isf_cyclic(i,:,j)'.*smoothvec'));
                    res(i,j).SKw = tmpSKw(indx)/max(tmpSKw)*max_intensity;
                    res(i,j).ISF0 = isf(i,1:decim_fctr:end,j);
                    res(i,j).max_ISF0 = max_isf(i,1,j);
                    res(i,j).setime = params.t_isf(1:decim_fctr:end);
                    res(i,j).mean.Preal = real(res(i,j).ISF0);
                    res(i,j).mean.Pimag = zeros(size(res(i,j).ISF0));
                    res(i,j).Energ_meV = Energ_meV;
                    res(i,j).gammaTrue = [];
                    res(i,j).E0 = [];
                    res(i,j).KxKy = params.dK(i,:,j);
                    % Determine azimuth
                    tmp = res(i,j).KxKy(1)/res(i,j).KxKy(2);
                    if abs(tmp) > 28, res(i,j).alpha=110; % within +/- 2deg
                    elseif abs(tmp-sqrt(3))<0.13, res(i,j).alpha=112; % within +/- 2deg
                    elseif isnan(tmp)
                        if j==1
                            res(i,j).alpha=110;
                        else
                            res(i,j).alpha=112;
                        end
                    else
                        res(i,j).alpha=112;
                    end
                    res(i,j).dK = sqrt(params.dK(i,1,j).*params.dK(i,1,j)+params.dK(i,2,j).*params.dK(i,2,j));
                    res(i,j).temperature = params.T;
                end
            end
            
            res = reshape(res,1,[]);
            
            % Sort by dK, is is more comfortable
            [tmp_dK,indx]=sort(abs([res.dK]));
            res = res(indx);
            
        end
        
        function res = fit_dyFiles_freq_domain(res,buildFitFunction,resCase,plotMode,refRes,axisize,fitrange)
            
            % ref Res is a res structure that provides the exponential
            % decay as fitted in the time domain. It can then be used to
            % constrain the fitting in SKw for t0+he jump part
            
            load_chess_parameters
            
            % Sort by dK, is is more comfortable (for the res file its done
            % in the preparation function
            if exist('refRes', 'var') && ~isempty(refRes)
                [tmp_dK,indx]=sort(abs([refRes.dK]));
                refRes = refRes(indx);
            end
            
            %% Fit
            isGraphicOn = usejava('desktop');
            if isGraphicOn, h = waitbar(0,'Please wait...'); end
            for j=1:length(res)
                
                isExpData= isfield(res(j),'endStatus');
                
                refRes_ = [];
                if exist('refRes', 'var') && ~isempty(refRes),refRes_ = refRes(j); end
                
                
                Energ_meV = res(j).Energ_meV; if size(Energ_meV,1) < size(Energ_meV,2), Energ_meV=Energ_meV'; end
                SKw_orig = real(res(j).SKw); if size(SKw_orig,1) < size(SKw_orig,2), SKw_orig=SKw_orig'; end
                
                
                % if experimental data, Apply detailed balance
                if ~isempty(res(j).E0)
                    indx = Energ_meV>0;
                    SKw_orig(indx) = SKw_orig(indx).*exp(Energ_meV(indx)/(SE_kB*res(j).temperature/(SE_e/1000)));
                end
                
                % ====== Construct the fitting parameters (Boundaries, starting guess, etc.) ======
                
                % find the elastic bin, then remove it (constant of the polarization curve)
                SKw = SKw_orig;
                elastic_bin = find(abs(Energ_meV)<1e-6);
                if isempty(elastic_bin), elastic_bin=find(abs(Energ_meV)-min(abs(Energ_meV))<=0); end
                elastSKw = mean(SKw(elastic_bin));
                SKw(elastic_bin)=[]; Energ_meV(elastic_bin)=[];
                
                % remove neighbourhood of elastic bin, to help fitting the
                % foot (comment if not desired)
                if 0
                    near_elastic_bins = abs(Energ_meV)<0.18; % removing ueVs
                    near_elastSKw = SKw(near_elastic_bins); % not used
                    SKw(near_elastic_bins)=[]; Energ_meV(near_elastic_bins)=[];
                end
                
                % Normalize, we still keep elastSKw to be able to scale
                SKw = SKw/elastSKw;
                
                % Get fitting function etc.
                [res(j).fncName,fitFuncStr,fitParams] = buildFitFunction(resCase, res(j).E0,res(j).temperature,res(j).dK, res(j).gammaTrue, 44.4,Energ_meV,SKw, refRes_);
                
                diffBoundaries = fitParams(:,2)-fitParams(:,1); tmpIndx = diffBoundaries<2e-8;
                fitParams(tmpIndx,2) = fitParams(tmpIndx,1)+2e-8;
                
                opts = fitoptions('method','nonLinearLeastSquares','Robust', 'off', 'lower', fitParams(:,1), 'upper', fitParams(:,2), 'startpoint',fitParams(:,3) , ...
                    'maxiter',6000,'MaxFunEvals',6000,'algorithm','trust-region');
                ftype = fittype(fitFuncStr,'options',opts);
                
                % ====== Check for the best (start) cut-off time ======
                % Fit to different start times, and work out the best rmse:
                % Find cut-off dE. For simulation we want only the QHAS
                % part (and a little bit of edges, as there is also the
                % T-mode inellastic peak. For experimental data, look for a
                % range that is best for fitting (so ignore for example the
                % jacobian amplification of the noise at negative end of
                % the spectrum, or don't take too much of the positive
                % side, to avoid enhancing unimportant parts in the
                % fitting).
                
                % set limits - first two are boundaries for the
                % minCutOffVec, third is for maxCutOffVec
                experLimits = [-5 -4 60];
                simulLimits = [-20.1 -20 20];
                
                if ~isExpData, limits = simulLimits; lngthOfLimitVecs = 1; else limits = experLimits; lngthOfLimitVecs = 2; end
                
                minCutOffVec = find(limits(1) < Energ_meV & Energ_meV <= limits(2));
                if isempty(minCutOffVec), minCutOffVec=1; end
                minCutOffVec = minCutOffVec(1:floor(length(minCutOffVec)/lngthOfLimitVecs)+1:end);
                maxCutOffVec = find(Energ_meV <= limits(3),1,'last');
                
                rmse = nan(length(minCutOffVec),length(maxCutOffVec));
                for minIndx=1:length(minCutOffVec)
                    if length(minCutOffVec) == 1 && length(maxCutOffVec) == 1, rmse = 1; break; end
                    parfor maxIndx=1:length(maxCutOffVec)
                        if maxCutOffVec(maxIndx)-minCutOffVec(minIndx) < 20
                            continue;
                        end
                        indxVec = minCutOffVec(minIndx):maxCutOffVec(maxIndx);
                        tmp_Energ_meV = Energ_meV(indxVec); tmp_SKw = SKw(indxVec);
                        [~, gof, ~] = fit(tmp_Energ_meV,tmp_SKw,ftype);
                        rmse(minIndx,maxIndx) = gof.rmse;
                    end
                end
                
                % ====== fit on the best range ======
                [r,c] = find(rmse==min(rmse(:)));
                res(j).dErange4fitting = [minCutOffVec(r):maxCutOffVec(c)];
                
                % if a fit range is supplied, use this range as the fitting
                % range.
                if exist('fitrange','var')
                    if length(fitrange)==2 && numel(fitrange)==2
                        res(j).dErange4fitting = find(fitrange(1) < Energ_meV & Energ_meV <= fitrange(2));
                    end
                end
                
                [ffun, gof, ~] = fit(Energ_meV(res(j).dErange4fitting),SKw(res(j).dErange4fitting),ftype);
                res(j).ffun = ffun;
                res(j).gof = gof;
                res(j).ci = confint(ffun,0.66);
                res(j).yBest=ffun(Energ_meV(res(j).dErange4fitting));
                res(j).vBest = coeffvalues(ffun);
                
                if plotMode
                    tmpMod=mod(j,6);
                    if tmpMod-1 == 0
                        fig1 = figure;
                        if isGraphicOn, set(fig1, 'WindowStyle', 'docked'); end
                    end
                    if tmpMod==0, tmpMod=6; end
                    subplot(3,2,tmpMod)
                    plot(res(j).Energ_meV,real(SKw_orig)/elastSKw)
                    x_orig=Energ_meV(res(j).dErange4fitting);
                    x=linspace(x_orig(1),x_orig(end),length(x_orig)*6);
                    fit_dyFiles.plot_funcional(x,ffun,gof,res(j).fncName);
                    title([res(j).filename(end-11:end-4) ' dK=' sprintf('%0.3f',res(j).dK) ' ' char(197) '^{-1}' ' T=' sprintf('%d',res(j).temperature) ' K'])
                    if exist('axisize','var'),if length(axisize)==4 && numel(axisize)==4, axis(axisize); end; end
                end
                
                if isGraphicOn, waitbar(j/length(res)); else disp(['finished res(' num2str(j/length(res)) ')']); end
                if isExpData
                    if isfield(res(j).mean,'bkg') && isfield(res(j).mean,'C')
                        res(j).elastpeak=elastSKw*(mean(res(j).mean.C)-res(j).mean.bkg);
                    end
                end
            end
            if isGraphicOn, close(h); end
        end
        
        function res = fit_dyFiles_time_domain(res, varargin)
            
            handPeak_tStart = 0;
            takeMagnitude = 1; % 0 - take Preal, 1 - take sqrt(Preal^2+Pimag^2), 2 - take magnitude but approximate Pimag with linear fit.
            disp(['takeMagnitude was set to ' num2str(takeMagnitude)])
            
            args4Fit_Func = {};
            if nargin > 1
                Fit_func = varargin{1};
                if nargin > 2
                    args4Fit_Func = varargin(2:end);
                end
            else
                Fit_func = @Fit.nExponentDecayFit;
            end
            
            if nargout < 1, disp('no res output param was defined, EXITING'); return; end
            
            load_chess_parameters
            
            %% Fit
            h = waitbar(0,'Please wait...');
            for j=1:length(res)
                clear startTime vBest yBest ci rsquare
                setime = res(j).setime; Pimag = res(j).mean.Pimag; Preal = res(j).mean.Preal;
                if size(setime,1)>size(setime,2), setime=setime'; end
                if size(Pimag,1)>size(Pimag,2), Pimag=Pimag'; end
                if size(Preal,1)>size(Preal,2), Preal=Preal'; end
                
                % If both positive and negative time exists ==> mirror
                if setime(1)<0
                    indx1=setime<0; setime(indx1)=[];
                    Preal1=fliplr(Preal(indx1)); Preal(indx1)=[]; Preal(2:end)=(Preal(2:end)+Preal1)/2;
                    Pimag1=fliplr(Pimag(indx1)); Pimag(indx1)=[]; Pimag(2:end)=(Pimag(2:end)+Pimag1)/2;
                end
                
                % get a subset of the ISF, based on min/max values
                indx = Preal>0.01;
                ind_indx = find(~indx,1,'first'); if isempty(ind_indx), ind_indx = length(indx); end
                indx = indx(1:ind_indx);
                if sum(indx) < 2000
                    indx = 1:length(Preal);
                end
                Preal = Preal(indx); Pimag = Pimag(indx); setime = setime(indx);
                
                isExpData= isfield(res(j),'endStatus');
                
                tStartRange=[2,min(0.5*max(setime),6)]; tIndxLength=5;
                %tStartRange=[30,min(0.5*max(setime),40)]; tIndxLength=2;
                if tStartRange(1)>tStartRange(2), tStartRange(2)=max(setime); end
                tStartIndxVec = find(setime > tStartRange(1) & setime < tStartRange(2));
                tStartIndxVec = tStartIndxVec(1:ceil(length(tStartIndxVec)/tIndxLength):length(tStartIndxVec));
                tStartIndxVec = unique(tStartIndxVec);
                if isempty(tStartIndxVec) || length(setime)-tStartIndxVec(end) < 20
                    tStartIndxVec = 1;
                end
                
                % Check for the best (start) cut-off time. Fit to different
                % start times, and work out the best rsquare.
                rsquare = ones(length(tStartIndxVec),1);
                for i=1:length(tStartIndxVec)
                    tStartIndx = tStartIndxVec(i);
                    
                    if takeMagnitude == 0
                        realP = Preal(tStartIndx:end);
                    elseif takeMagnitude == 1
                        realP = sqrt(Preal(tStartIndx:end).^2+Pimag(tStartIndx:end).^2);
                    elseif takeMagnitude == 2
                        [~,yBestImag,~,~] = Fit.linearFit(setime(tStartIndx:end),Pimag(tStartIndx:end));
                        realP = sqrt(Preal(tStartIndx:end).^2+yBestImag.^2);
                    end
                    [~,~,~,rsquare(i)] = Fit_func(setime(tStartIndx:end),realP,args4Fit_Func{:});
                end
                
                if handPeak_tStart
                    h1 = figure();
                    subplot(2,1,1); xlabel('tStart'); ylabel('rsquare'); plot(setime(tStartIndxVec),rsquare,'o'); set(gca, 'YScale', 'log')
                    subplot(2,1,2); xlabel('t_{SE}'); ylabel('ISF'); plot(setime,Preal,'o'); set(gca, 'YScale', 'log')
                    chosenPoints = CustomFuncs.peakPointsFromFigure('dataX',setime(tStartIndxVec),'dataY',rsquare);
                    close(h1)
                    indx = find(chosenPoints,1,'first');
                    if isempty(indx), indx = 1; end
                else
                    [~,indx] = max(rsquare);
                end
                res(j).startTime = setime(tStartIndxVec(indx));
                res(j).endTime = setime(end);
                setimeIndx = setime > res(j).startTime;
                
                res(j).takeMagnitude = takeMagnitude;
                if takeMagnitude
                    [res(j).vBestImag,res(j).yBestImag,res(j).ciImag,res(j).rsquareImag] = Fit.linearFit(setime(setimeIndx),Pimag(setimeIndx));
                    res(j).realP = sqrt(Preal(setimeIndx).^2+res(j).yBestImag.^2);
                else
                    res(j).realP = Preal(setimeIndx);
                end
                [res(j).vBest,res(j).yBest,res(j).ci,res(j).rsquare,res(j).ffun] = Fit_func(setime(setimeIndx),res(j).realP,args4Fit_Func{:});
                waitbar(j/length(res))
            end
            close(h)
            %save('tmp_res.mat','res')
        end
        
        function plot_single_SKw(res)
            % plot the ISF and fit for each measurement (or each combined measurement files)
            for j=1:length(res)
                if mod(j-1,6) == 0
                    h=figure; hold on;
                    if  usejava('desktop'), set(h, 'WindowStyle', 'docked'); end
                end
                figure(h); subplot(3,2,mod(j-1,6)+1); hold on;
                SKw_real = real(res(j).SKw); I_elastic = SKw_real(find(abs(res(j).Energ_meV)<1e-6,1));
                plot(res(j).Energ_meV,SKw_real/I_elastic)
                if isfield(res(j),'fncName'),
                    fit_dyFiles.plot_funcional(res(j).Energ_meV(res(j).dErange4fitting),res(j).ffun,res(j).gof,res(j).fncName);
                end
                if ~isfield(res(j),'alpha'), res(j).alpha=[];end
                title([res(j).filename(end-11:end-4) ' dK=' sprintf('%0.3f',res(j).dK)  '[' char(197) '^{-1}]' ...
                    ' T=' num2str(res(j).temperature) ' alpha=' num2str(res(j).alpha)])
            end
        end
        
        function plot_single_ISF(res,varargin)
            
            % parse varargin
            prsdArgs = inputParser;   % Create instance of inputParser class.
            prsdArgs.addParameter('logPlot',1, @isnumeric);
            prsdArgs.addParameter('plotReal',1, @isnumeric);
            prsdArgs.addParameter('plotImag',1, @isnumeric);
            prsdArgs.parse(varargin{:});
            logPlot = prsdArgs.Results.logPlot;
            plotReal = prsdArgs.Results.plotReal;
            plotImag = prsdArgs.Results.plotImag;
            
            % plot the ISF and fit for each measurement (or each combined measurement files)
            for j=1:length(res)
                if mod(j-1,9) == 0, h=figure; hold on; end
                figure(h); subplot(3,3,mod(j-1,9)+1)
                
                setimeIndx = res(j).setime > res(j).startTime & res(j).setime <= res(j).endTime;
                setime_indxed = res(j).setime(setimeIndx);
                setime = res(j).setime;
                
                hold on
                if plotReal, plot(setime,res(j).mean.Preal,'go'); end
                if plotImag, plot(setime,res(j).mean.Pimag,'ro'); end
                plot(setime_indxed,res(j).realP,'b.')
                plot(setime_indxed,res(j).ffun(setime_indxed),'k-')
                
                if isfield(res,'yBestImag')
                    plot(setime_indxed,res(j).yBestImag,'k-');
                end
                
                % plot individual exponentials
                numOfExp = (length(coeffnames(res(j).ffun))-1)/2;
                if numOfExp > 1
                    e_ =  res(j).vBest(end);
                    for i=1:2:numOfExp*2
                        a_(i) =  res(j).vBest(i);
                        b_(i)  = res(j).vBest(i+1);
                        plot(setime,a_(i) .* exp(-b_(i)*setime)+e_,'.');
                    end
                end
                
                if isfield(res,'alpha'), alpha=num2str(res(j).alpha);
                elseif isfield(res,'endStatus'), alpha=num2str(res(j).endStatus.alpha);
                end
                title([res(j).filename(end-11:end-4) ' dK = ' num2str(res(j).dK) ' T=' num2str(res(j).temperature) ' alpha=' alpha])
                if logPlot, set(gca, 'XScale', 'log'); end
                
            end
        end
        
        function plot_grouped_data_SKw_Dk(res,cellOfDEstr)
            for i=1:length(res)
                coeff_str = coeffnames(res(i).ffun);
                coeff_val = coeffvalues(res(i).ffun);
                for j=1:length(cellOfDEstr)
                    tmp = strfind(coeff_str,cellOfDEstr(j));
                    for k=1:length(tmp), if tmp{k}==1, break; end; end
                    if k==length(tmp), continue; end
                    dE(i,j) = coeff_val(k);
                    dKarray(i,j) = de2dk(res(i).E0,dE(i,j),res(i).gammaTrue,44.4);
                    %if coeff_val(k-1) < 0.02 || 0, dE(i,j) = NaN; dKarray(i,j) = NaN; end
                end
            end
            %             dE = dE';
            %             dKarray = dKarray';
%             figure;
            plot(abs(dKarray),abs(dE),'*','color','b');
            xlabel(join(['\DeltaK/',char(197),'^{-1}']));
        end
        
        function plot_grouped_data_alphaDk(res,varargin)
            % This function will collect all the fitting values form a
            % 'res' structure-array, and will plot the requested ones.
            % For each res(j), we expect to have a fitting string that can
            % be extracted from res(j).ffun, and a series of mode-names
            % which are held in res(j).fncName. For example,
            % 'a1*exp(-a2*x)+a3+b1*exp(-x^2/b2)' will be accompanied by
            % fncName={dcy_SOMENAME cg_SOMEOTHERNAME}. The refix will say
            % what type of functional is expected, and the SOMENAME will
            % say what is it used for (usually).
            % Known prefix: dcy - exponential decay, cl_/cg_/ncl/ncg -
            % centered (or non-centered) lorentzian or gaussian, sp1 - a
            % functional for T-modes.
            % The function doesn't assume that all res(j) have the same
            % fitting string.
            
            
            %% parse varargin
            prsdArgs = inputParser;   % Create instance of inputParser class.
            prsdArgs.addParameter('domain', 'dE', @ischar);% Values can be 'dE' or 't'.
            prsdArgs.addParameter('functionals', 'cl cg', @ischar);% Values can be: cl, cg, ncl, ncg, dcy
            prsdArgs.addParameter('rsquareCutOff', 0.8, @isnumeric);% Values can be: cl, cg, ncl, ncg, dcy
            prsdArgs.parse(varargin{:});
            
            domain = prsdArgs.Results.domain;
            functional_names = strsplit(prsdArgs.Results.functionals);
            rsquareCutOff = prsdArgs.Results.rsquareCutOff;
            
            
            % remove cases with rsquare < rsquareCutOff
            if isfield(res,'gof')
                for j=1:length(res), rsquare(j)=res(j).gof.rsquare; end
            elseif isfield(res,'rsquare')
                rsquare = [res.rsquare];
            end
            
            res(rsquare < rsquareCutOff) = [];
            
            %% Set variables for general use
            load_chess_parameters;
            isGraphicOn = usejava('desktop');
            
            uniq_specAtten = [2,10,50:50:1000];
            
            for j=1:length(res),
                if ~isfield(res(j),'clean_spec')
                    clean_spec_counts = 4.4e-8 %input('Whats the clean specular?');
                    break;
                end
            end
            
            for j=1:length(res)
                %% Set some values in res structure - this should be done elsewhere, in the creation of 'res' structure
                % Get specular counts
                if ~isfield(res(j),'clean_spec')
                    res(j).clean_spec_counts = clean_spec_counts;
                end
                
                % ========   If its possible to calculate attenuation - why not? ========
                if isfield(res(j),'spec')
                    % Assign the specular counts to a field that can be used to
                    % form matrix based assignements
                    res(j).spec_counts = res(j).spec.counts;
                else
                    res(j).spec_counts = NaN;
                end
                
                % Calculate how much specular was attenua
                spec_rel_cnts = [res(j).clean_spec_counts]./[res(j).spec_counts];
                [tmp ind] = min(abs(spec_rel_cnts-uniq_specAtten));
                res(j).spec_atten = uniq_specAtten(ind);
                
                if strcmp(domain,'dE')
                    % ======== Calculate the SKw normalization factor ========
                    % This value should really be stored in the fitting
                    % procedure, so this is a temporary solution.
                    elastic_bin = find(abs(res(j).Energ_meV)<1e-6);
                    if isempty(elastic_bin), elastic_bin=find(abs(res(j).Energ_meV)-min(abs(res(j).Energ_meV))<=0); end
                    res(j).elastSKw = real(res(j).SKw(elastic_bin));
                end
            end
            
            temperature = [res.temperature];
            if isfield(res,'alpha'), alpha = [res.alpha];
            elseif isfield(res,'endStatus')
                endStatus_strct_array = [res.endStatus];
                for alpha_ind = 1:length(endStatus_strct_array)
                    if isempty(endStatus_strct_array(alpha_ind).alpha)
                        endStatus_strct_array(alpha_ind).alpha = endStatus_strct_array(alpha_ind-1).alpha;
                    end
                end
                alpha = [endStatus_strct_array.alpha];
            end
            uniqT = unique(temperature);
            uniqAlpha = unique(alpha);
            
            fig1 = figure;
            if isGraphicOn, set(fig1, 'WindowStyle', 'docked'); end
            if ~isGraphicOn, set(fig1, 'Visible', 'off'); end
            for alpha_indx=1:length(uniqAlpha)
                if isnan(uniqAlpha(alpha_indx)), continue; end % dK=zero, both azimuth
                if alpha_indx > 1 && 0
                    fig1 = figure();
                    if  isGraphicOn, set(fig1, 'WindowStyle', 'docked'); end
                    if ~isGraphicOn, set(fig1, 'Visible', 'off'); end
                end
                for T_indx=1:length(uniqT)
                    indx = temperature == uniqT(T_indx) & (alpha == uniqAlpha(alpha_indx));
                    res_ = res(indx);
                    vBest=[res_.vBest];
                    ci = [res_.ci];
                    dK=[res_.dK];
                    num1 = length(vBest)/length(res_);
                    num2 = length(res_)-1;
                    coeffName=coeffnames(res_(1).ffun);
                    
                    switch domain
                        case 'dE'
                            elastSKw=[res_.elastSKw];
                            if size(elastSKw,1) > 1, elastSKw = mean(elastSKw,1); end; %if (by a mistake) more than one value was included for each dK
                            
                            a1_i = find(strcmp(coeffnames(res(1).ffun),'a1'),1);
                            a2_i = find(strcmp(coeffnames(res(1).ffun),'a2'),1);
                            a3_i = find(strcmp(coeffnames(res(1).ffun),'a3'),1);
                            a4_i = find(strcmp(coeffnames(res(1).ffun),'a4'),1);
                            b1_i = find(strcmp(coeffnames(res(1).ffun),'b1'),1);
                            b2_i = find(strcmp(coeffnames(res(1).ffun),'b2'),1);
                            b3_i = find(strcmp(coeffnames(res(1).ffun),'b3'),1);
                            b4_i = find(strcmp(coeffnames(res(1).ffun),'b4'),1);
                            
                            clear a1 a1_ci a2 a2_ci a3 a3_ci a4 a4_ci b1 b1_ci b2 b2_ci b3 b3_ci b4 b4_ci
                            for i=1:length(res_)
                                if ~isempty(a1_i), a1_ci(:,i)=res_(i).ci(:,a1_i); a1(i)=res_(i).vBest(:,a1_i); end
                                if ~isempty(a2_i), a2_ci(:,i)=res_(i).ci(:,a2_i); a2(i)=res_(i).vBest(:,a2_i); end
                                if ~isempty(a3_i), a3_ci(:,i)=res_(i).ci(:,a3_i); a3(i)=res_(i).vBest(:,a3_i); end
                                if ~isempty(a4_i), a4_ci(:,i)=res_(i).ci(:,a4_i); a4(i)=res_(i).vBest(:,a4_i); end
                                if ~isempty(b1_i), b1_ci(:,i)=res_(i).ci(:,b1_i); b1(i)=res_(i).vBest(:,b1_i); end
                                if ~isempty(b2_i), b2_ci(:,i)=res_(i).ci(:,b2_i); b2(i)=res_(i).vBest(:,b2_i); end
                                if ~isempty(b3_i), b3_ci(:,i)=res_(i).ci(:,b3_i); b3(i)=res_(i).vBest(:,b3_i); end
                                if ~isempty(b4_i), b4_ci(:,i)=res_(i).ci(:,b4_i); b4(i)=res_(i).vBest(:,b4_i); end
                            end
                            
                            if ~isempty(a1_i),
                                lorentzian_1_int_ci = a1_ci.*repmat(elastSKw,2,1); lorentzian_1_fwhm_ci = a2_ci;
                                lorentzian_1_int = a1.*elastSKw; lorentzian_1_fwhm = a2;
                                exp_dcy_1_ci = a2_ci/1000*SE_e*pi/SE_h/1e9;
                                exp_dcy_1 = a2/1000*SE_e*pi/SE_h/1e9;
                            end
                            
                            if ~isempty(a3_i)
                                lorentzian_2_int_ci = a3_ci.*repmat(elastSKw,2,1); lorentzian_2_fwhm_ci = a4_ci;
                                lorentzian_1_int = a3.*elastSKw; lorentzian_1_fwhm = a4;
                                exp_dcy_2_ci = a4_ci/1000*SE_e*pi/SE_h/1e9;
                                exp_dcy_2 = a4/1000*SE_e*pi/SE_h/1e9;
                            end
                            
                            if ~isempty(b1_i)
                                gaussian_1_int_ci = b1_ci.*b2_ci*sqrt(pi).*repmat(elastSKw,2,1); gaussian_1_fwhm_ci = b2_ci*2*sqrt(log(2));
                                gaussian_1_int = b1.*b2*sqrt(pi).*elastSKw; gaussian_1_fwhm = b2*2*sqrt(log(2));
                            end
                            
                            if ~isempty(b3_i)
                                gaussian_2_int_ci = b3_ci.*b4_ci*sqrt(pi).*repmat(elastSKw,2,1); gaussian_2_fwhm_ci = b4_ci*2*sqrt(log(2));
                                gaussian_2_int = b3.*b4*sqrt(pi).*elastSKw; gaussian_2_fwhm = b4*2*sqrt(log(2));
                            end
                            
                            if exist('lorentzian_1_int') && exist('gaussian_1_int')
                                rel_lorentz_1 = lorentzian_1_int./(lorentzian_1_int+gaussian_1_int);
                            end
                            if exist('lorentzian_2_int') && exist('gaussian_2_int')
                                rel_lorentz_2 = lorentzian_2_int./(lorentzian_2_int+gaussian_2_int);
                            end
                            
                            subplot(3,1,1); hold on;
                            if exist('exp_dcy_1'), errorbar(abs(dK),exp_dcy_1,diff(exp_dcy_1_ci)/2,'*'); end
                            if exist('exp_dcy_2'), errorbar(abs(dK),exp_dcy_2,diff(exp_dcy_2_ci)/2,'*'); end
                            xlabel(['\Delta K ' '[' char(197) '^{-1}]'],'FontSize',14)
                            ylabel('\alpha [GHz]','FontSize',14)
                            title(['Filenames: ' res_(1).filename ' - ' res_(end).filename '  ... T=' num2str(uniqT(T_indx)) '[K] \alpha=' num2str(uniqAlpha(alpha_indx)) ' I_0/I=']);
                            
                            subplot(3,1,2); hold on;
                            if exist('gaussian_1_fwhm'), errorbar(abs(dK),gaussian_1_fwhm,diff(gaussian_1_fwhm_ci)/2,'o'); end
                            if exist('gaussian_2_fwhm'), errorbar(abs(dK),gaussian_2_fwhm,diff(gaussian_2_fwhm_ci)/2,'o'); end
                            xlabel(['\Delta K ' '[' char(197) '^{-1}]'],'FontSize',14)
                            ylabel('FWHM [meV]','FontSize',14)
                            title(['T=' '[K] \alpha=' ' I_0/I=']);
                            
                            subplot(3,1,3); hold on;
                            if exist('rel_lorentz_1'), errorbar(abs(dK),rel_lorentz_1,zeros(size(rel_lorentz_1)),'b*'); end
                            if exist('rel_lorentz_2'), errorbar(abs(dK),rel_lorentz_2,zeros(size(rel_lorentz_2)),'b*'); end
                            xlabel(['\Delta K ' '[' char(197) '^{-1}]'],'FontSize',14)
                            ylabel({'Jump Integral','-------------------','(Jump+Ballistic integrals)'},'FontSize',14)
                            title(['T=' '[K] \alpha=' ' I_0/I=']);
                            
                        case 't'
                            
                            P0_ci=[]; P0=[];
                            exp_dcy_ci=[]; exp_dcy=[];
                            cnst_ofst_ci=[]; cnst_ofst=[];
                            
                            argNames = argnames(res_(1).ffun);
                            nExponent = 0;
                            while argNames{nExponent+1}(1)=='a'
                                nExponent = nExponent + 1;
                            end
                            
                            isCnstOfst = 0;
                            if find([argNames{:}]=='e',1,'first'), isCnstOfst = 1; end
                            
                            tmp_indx=1;
                            
                            for exp_indx=1:nExponent % num of exponents in slow decay
                                P0_ci(:,:,exp_indx)=ci(:,tmp_indx:num1:num2*num1+tmp_indx);
                                P0(:,:,exp_indx)=vBest(:,tmp_indx:num1:num2*num1+tmp_indx);
                                tmp_indx = tmp_indx + 1;
                                exp_dcy_ci(:,:,exp_indx)=ci(:,tmp_indx:num1:num2*num1+tmp_indx);
                                exp_dcy(:,:,exp_indx)=vBest(:,tmp_indx:num1:num2*num1+tmp_indx);
                                tmp_indx = tmp_indx + 1;
                            end
                            
                            if isCnstOfst
                                cnst_ofst_ci(:,:,exp_indx)=ci(:,tmp_indx:num1:num2*num1+tmp_indx);
                                cnst_ofst(:,:,exp_indx)=vBest(:,tmp_indx:num1:num2*num1+tmp_indx);
                            end
                            
                            exp_dcy_ci = exp_dcy_ci * 1000; % GHz
                            exp_dcy = exp_dcy * 1000; % GHz
                            
                            Legend{alpha_indx} = ['\alpha=' num2str(uniqAlpha(alpha_indx))];
                            subplot(3,1,1); hold on;
                            for tmp_indx=1:nExponent
                                errorbar(abs(dK),squeeze(exp_dcy(:,:,tmp_indx)),diff(squeeze(exp_dcy_ci(:,:,tmp_indx)))/2,'*');
                            end
                            xlabel(['\Delta K ' '[' char(197) '^{-1}]'],'FontSize',14)
                            ylabel('\alpha [GHz]','FontSize',14)
                            title(['T=' num2str(uniqT(T_indx)) '[K] I_0/I=']);
                            legend(Legend)
                            
                            subplot(3,1,2); hold on;
                            for tmp_indx=1:nExponent
                                errorbar(abs(dK),P0(:,:,tmp_indx),diff(P0_ci(:,:,tmp_indx))/2,'*');
                            end
                            xlabel(['\Delta K ' '[' char(197) '^{-1}]'],'FontSize',14)
                            ylabel('P0','FontSize',14)
                            legend(Legend)
                            
                            if isCnstOfst
                                subplot(3,1,3); hold on;
                                for tmp_indx=1:nExponent
                                    errorbar(abs(dK),cnst_ofst(:,:,tmp_indx),diff(cnst_ofst_ci(:,:,tmp_indx))/2,'*');
                                end
                                xlabel(['\Delta K ' '[' char(197) '^{-1}]'],'FontSize',14)
                                ylabel('Constant','FontSize',14)
                                legend(Legend)
                            end
                            
                        otherwise
                            error('no such domain')
                    end
                end
                saveas(fig1, ['Untitled_tmp_' num2str(fig1.Number)],'fig')
            end
        end
        
        function plot_grouped_data_arrhenius(res,varargin)
            % This function will collect all the fitting values form a
            % 'res' structure-array, and will plot the requested ones.
            % For each res(j), we expect to have a fitting string that can
            % be extracted from res(j).ffun, and a series of mode-names
            % which are held in res(j).fncName. For example,
            % 'a1*exp(-a2*x)+a3+b1*exp(-x^2/b2)' will be accompanied by
            % fncName={dcy_SOMENAME cg_SOMEOTHERNAME}. The refix will say
            % what type of functional is expected, and the SOMENAME will
            % say what is it used for (usually).
            % Known prefix: dcy - exponential decay, cl_/cg_/ncl/ncg -
            % centered (or non-centered) lorentzian or gaussian, sp1 - a
            % functional for T-modes.
            % The function doesn't assume that all res(j) have the same
            % fitting string.
            
            
            %% parse varargin
            prsdArgs = inputParser;   % Create instance of inputParser class.
            prsdArgs.addParameter('domain', 'dE', @ischar);% Values can be 'dE' or 't'.
            prsdArgs.addParameter('functionals', 'cl cg', @ischar);% Values can be: cl, cg, ncl, ncg, dcy
            prsdArgs.addParameter('tempRange', [], @isnumeric);% [min max]
            prsdArgs.parse(varargin{:});
            
            domain = prsdArgs.Results.domain;
            functional_names = strsplit(prsdArgs.Results.functionals);
            tempRange = prsdArgs.Results.tempRange;
            
            
            %% Set variables for general use
            load_chess_parameters;
            isGraphicOn = usejava('desktop');
            
            uniq_specAtten = [2,10,50:50:1000];
            
            for j=1:length(res),
                if ~isfield(res(j),'clean_spec')
                    clean_spec_counts = 4.4e-8 %input('Whats the clean specular?');
                    break;
                end
            end
            
            for j=1:length(res)
                %% Set some values in res structure - this should be done elsewhere, in the creation of 'res' structure
                % Get specular counts
                if ~isfield(res(j),'clean_spec')
                    res(j).clean_spec_counts = clean_spec_counts;
                end
                
                % ========   If its possible to calculate attenuation - why not? ========
                if isfield(res(j),'spec')
                    % Assign the specular counts to a field that can be used to
                    % form matrix based assignements
                    res(j).spec_counts = res(j).spec.counts;
                else
                    res(j).spec_counts = NaN;
                end
                
                % Calculate how much specular was attenua
                spec_rel_cnts = [res(j).clean_spec_counts]./[res(j).spec_counts];
                [tmp ind] = min(abs(spec_rel_cnts-uniq_specAtten));
                res(j).spec_atten = uniq_specAtten(ind);
                
                if strcmp(domain,'dE')
                    % ======== Calculate the SKw normalization factor ========
                    % This value should really be stored in the fitting
                    % procedure, so this is a temporary solution.
                    elastic_bin = find(abs(res(j).Energ_meV)<1e-6);
                    res(j).elastSKw = real(res(j).SKw(elastic_bin));
                end
                
            end
            
            dK = [res.dK];
            if isfield(res,'alpha'), alpha = [res.alpha];
            elseif isfield(res,'endStatus')
                endStatus_strct_array = [res.endStatus];
                for alpha_ind = 1:length(endStatus_strct_array)
                    if isempty(endStatus_strct_array(alpha_ind).alpha)
                        if alpha_ind == 1
                            endStatus_strct_array(alpha_ind).alpha = endStatus_strct_array(alpha_ind+1).alpha;
                        else
                            endStatus_strct_array(alpha_ind).alpha = endStatus_strct_array(alpha_ind-1).alpha;
                        end
                    end
                end
                alpha = [endStatus_strct_array.alpha];
            end
            uniq_dK = uniquetol(dK,0.02);
            uniqAlpha = uniquetol(alpha,1);
            
            fig1 = figure;
            if  isGraphicOn, set(fig1, 'WindowStyle', 'docked'); end
            for alpha_indx=1:length(uniqAlpha)
                if isnan(uniqAlpha(alpha_indx)), continue; end % dK=zero, both azimuth
                if alpha_indx > 1 && 0
                    fig1 = figure;
                    if  isGraphicOn, set(fig1, 'WindowStyle', 'docked');end
                    subplot(3,1,1);
                    title(['Filenames: ' res_(1).filename ' - ' res_(end).filename '  ... dK=' num2str(uniq_dK(dK_indx)) '[' char(197) '^{-1}]' ' \alpha=' num2str(uniqAlpha(alpha_indx)) ' I_0/I=']);
                    subplot(3,1,2);
                    title(['Filenames: ' res_(1).filename ' - ' res_(end).filename '  ... dK=' num2str(uniq_dK(dK_indx)) '[' char(197) '^{-1}]' ' \alpha=' num2str(uniqAlpha(alpha_indx)) ' I_0/I=']);
                    
                end
                
                oneOverT_legend = {}; lgnd={}; markers = '*sod<^*sod<^*sod<^'; colors = 'brmkbrmk';
                for dK_indx=1:length(uniq_dK)
                    indx = ismembertol(dK,uniq_dK(dK_indx),2.1e-2) & (ismembertol(alpha,uniqAlpha(alpha_indx),1.1) | isnan(alpha));
                    res_ = res(indx);
                    if exist('tempRange','var') && ~isempty(tempRange), res_ = res_( tempRange(1) <= [res_.temperature] & [res_.temperature] <= tempRange(2) ); end
                    [temperature, ind_by_temp]=sort([res_.temperature]);
                    res_ = res_(ind_by_temp);
                    vBest=[res_.vBest];
                    ci = [res_.ci];
                    
                    
                    num1 = length(vBest)/length(res_);
                    num2 = length(res_)-1;
                    
                    clear a_ci b_ci e_ci a1_ci a2_ci b1_ci b2_ci gaussian_int_ci gaussian_fwhm_ci lorentzian_int_ci lorentzian_fwhm_ci a_ a1 a2 b_ b1 b2 e_
                    
                    oneOverT_legend(dK_indx) = {['dK=' num2str(uniq_dK(dK_indx))]};
                    if strcmp(domain,'dE')
                        elastSKw=[res_.elastSKw];
                        for i=1:length(res_), a1_ci(:,i)=res_(i).ci(:,1);  a2_ci(:,i)=res_(i).ci(:,2); b1_ci(:,i)=res_(i).ci(:,3); b2_ci(:,i)=res_(i).ci(:,4); end
                        gaussian_int_ci = b1_ci.*b2_ci*sqrt(pi).*repmat(elastSKw,2,1); gaussian_fwhm_ci = b2_ci*2*sqrt(log(2)); lorentzian_int_ci = a1_ci.*repmat(elastSKw,2,1); lorentzian_fwhm_ci = a2_ci;
                        exp_dcy_ci = a2_ci/1000*SE_e*pi/SE_h/1e9;
                        
                        for i=1:length(res_), a1(i)=res_(i).vBest(:,1);  a2(i)=res_(i).vBest(:,2); b1(i)=res_(i).vBest(:,3); b2(i)=res_(i).vBest(:,4); end
                        gaussian_int = b1.*b2*sqrt(pi).*elastSKw; gaussian_fwhm = b2*2*sqrt(log(2)); lorentzian_int = a1.*elastSKw; lorentzian_fwhm = a2;
                        exp_dcy = a2/1000*SE_e*pi/SE_h/1e9;
                        
                        rel_lorentz = lorentzian_int./(lorentzian_int+gaussian_int);
                        
                        subplot(2,2,1); hold on;
                        errorbar(1./temperature,log(exp_dcy),diff(log(exp_dcy_ci))/2,'*-');
                        xlabel('Temperature [1/K]','FontSize',14)
                        ylabel('ln(\alpha [GHz])','FontSize',14)
                        legend(oneOverT_legend)
                        
                        subplot(2,2,2); hold on;
                        errorbar(1./temperature,rel_lorentz,zeros(size(rel_lorentz)),'b*');
                        xlabel('Temperature [1/K]','FontSize',14)
                        ylabel({'Jump Integral','-------------------','(Jump+Ballistic integrals)'},'FontSize',14)
                        legend(oneOverT_legend)
                        
                    elseif strcmp(domain,'t')
                        numOfExp = (length(coeffnames(res(1).ffun))-1)/2;
                        
                        for i=1:numOfExp
                            indStrt = (i-1)*2+1;
                            indVec = indStrt:num1:num2*num1+indStrt;
                            a_ci(:,:,i) = ci(:,indVec);
                            a_(:,i) =  vBest(:,indVec);
                            indStrt = (i-1)*2+2;
                            indVec = indStrt:num1:num2*num1+indStrt;
                            b_ci(:,:,i) = ci(:,indVec);
                            b_(:,i)  = vBest(:,indVec);
                        end
                        indStrt = numOfExp*2+1;
                        e_ci(:,:,i) = ci(:,indStrt:num1:num2*num1+indStrt);
                        e_(:,i) =  vBest(:,indStrt:num1:num2*num1+indStrt);
                        
                        b_ci = b_ci * 1000; % GHz
                        b_ = b_ * 1000; % GHz
                        
                        subplot(2,2,1); hold on;
                        for i=1:numOfExp
                            errorbar(1./temperature,log(b_(:,i)),diff(log(b_ci(:,:,i)))/2,[markers(dK_indx) colors(dK_indx)]);
                            lgnd = {lgnd{:},oneOverT_legend{dK_indx}};
                        end
                        xlabel('Temperature [1/K]','FontSize',14)
                        ylabel('ln(\alpha [GHz])','FontSize',14)
                        legend(lgnd)
                        
                    else
                        error('no such domain')
                    end
                    
                    % fit to linear
                    for i=1:numOfExp
                        invT = 1./temperature';
                        log_b = log(b_(:,i));
                        includeIndx = CustomFuncs.manuallyCleanNoisyDataPoints(invT,log_b);
                        ftype = fittype('A-B*x');
                        [ffun, gof, ~] = fit(invT(includeIndx),log_b(includeIndx),ftype);
                        A(dK_indx,i) = exp(ffun.A);
                        dEb(dK_indx,i) = ffun.B * SE_kB / SE_e * 1000; % meV
                        ci = confint(ffun,0.66);
                        A_ci(:,dK_indx,i) = ci(:,1);
                        dEb_ci(:,dK_indx,i) = ci(:,2) * SE_kB / SE_e * 1000; % meV
                        
                        hold on; plot(1./temperature,ffun(1./temperature))
                        
                    end
                end
                
                subplot(2,2,3); hold on;
                for i=1:numOfExp
                    A_ci_tmp = A_ci(:,:,i);
                    errorbar(uniq_dK,A(:,i),diff(A_ci_tmp)/2,'*');
                end
                xlabel('dK','FontSize',14)
                ylabel('Pre-Exp','FontSize',14)
                
                subplot(2,2,4); hold on;
                for i=1:numOfExp
                    dEb_ci_tmp = dEb_ci(:,:,i);
                    errorbar(uniq_dK,dEb(:,i),diff(dEb_ci_tmp)/2,'*');
                end
                xlabel('dK','FontSize',14)
                ylabel('Eb [meV]','FontSize',14)
            end
        end
        
        function plot_seperated_data(res)
            h_alpha_dK=zeros(2,1);
            for j=1:length(res)
                indxOfFigure=find(res(j).temperature == h_alpha_dK(1,:), 1);
                if isempty(indxOfFigure)
                    h_alpha_dK(1,size(h_alpha_dK,2)+1)=res(j).temperature;
                    h_alpha_dK(2,size(h_alpha_dK,2)) = figure;
                    indxOfFigure=size(h_alpha_dK,2);
                end
                
                figure(h_alpha_dK(2,indxOfFigure)); title(['Temperature:' num2str(h_alpha_dK(1,indxOfFigure))]); hold on;
                
                subplot(3,1,1); hold on;
                hnew = errorbar(res(j).dK+j*0,sum(res(j).ci(:,2))/2*res(j).factor2,diff(res(j).ci(:,2))/2*res(j).factor2);
                xlabel(['\Delta K ' '[' char(197) '^{-1}]'],'FontSize',14)
                %xlabel(['Number of measurement'],'FontSize',14)
                ylabel('\alpha [psec^{-1}]','FontSize',14)
                
                subplot(3,1,2); hold on;
                hnew = errorbar(res(j).dK+j*0,sum(res(j).ci(:,1))/2,diff(res(j).ci(:,1))/2);
                xlabel(['\Delta K ' '[' char(197) '^{-1}]'],'FontSize',14)
                %xlabel(['Number of measurement'],'FontSize',14)
                ylabel('P_0','FontSize',14)
                
                subplot(3,1,3); hold on;
                hnew = errorbar(res(j).dK+j*0,sum(res(j).ci(:,3))/2,diff(res(j).ci(:,3))/2);
                xlabel(['\Delta K ' '[' char(197) '^{-1}]'],'FontSize',14)
                %xlabel(['Number of measurement'],'FontSize',14)
                ylabel('Constant','FontSize',14)
                
                legend('show')
                if strcmp(version,'8.4.0.150421 (R2014b)')
                    legObj = legend;
                    if strcmp(legObj.String,'data1')
                        legend_str = {[res(j).filename ' ' num2str(res(j).temperature) ' dK=' num2str(res(j).dK)]};
                    else
                        legend_str = {legObj.String,[res(j).filename ' ' num2str(res(j).temperature) ' dK=' num2str(res(j).dK)]};
                    end
                    legend(legend_str{:});
                else
                    % Get object handles
                    [LEGH,OBJH,OUTH,OUTM] = legend;
                    % Add object with new handle and new legend string to legend
                    legend([OUTH;hnew],OUTM{:},[res(j).filename ' ' num2str(res(j).temperature) ' dK=' num2str(res(j).dK)])
                end
            end
        end
        
        function [ffun,gof] = fit_Chadely_Elliot(dK,exp_dcy,azimuth)
            
            if ~size(dK,1)<size(dK,2), dK=dK'; end
            if ~size(exp_dcy,1)<size(exp_dcy,2), exp_dcy=exp_dcy'; end
            
            % ==== Fit to Chadley Elliot ====
            dK4fit = dK(dK<3); exp_dcy4fit = exp_dcy(dK<3);
            ftype = fittype(['BProb_baravis(2.5546,a1,x,[],' num2str(azimuth) ')'],'coef',{'a1'});
            [ffun_tmp, ~, ~] = fit(dK4fit,exp_dcy4fit,ftype);
            indx4fit = abs(ffun_tmp(dK4fit)-exp_dcy4fit) < ffun_tmp(dK4fit)*2;
            [ffun, gof, ~] = fit(dK4fit(indx4fit),exp_dcy4fit(indx4fit),ftype);
        end
        
    end
end
