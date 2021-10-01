classdef postprocess_dyfiles

        %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(Constant)
        projFold='Cu111_4';
        %projFold='/home/na364/Downloads/H2OORu/Temperature/';
    end
    
    methods(Static)
        
        function postprocess_dyfile(varargin)
            
            % parse varargin
            prsdArgs = inputParser;   % Create instance of inputParser class.
            prsdArgs.addParameter('filenames', {}, @iscell);
            prsdArgs.addParameter('initStr', 'dy0', @ischar);
            prsdArgs.addParameter('NumVec', [], @isnumeric);
            prsdArgs.addParameter('reProcessing', 0, @isnumeric);
            prsdArgs.addParameter('fixFalsePositive', 0, @isnumeric);
            prsdArgs.addParameter('keepIlength', 0, @isnumeric);
            prsdArgs.addParameter('intpI', [], @isnumeric);
            prsdArgs.parse(varargin{:});
            
            filenames = prsdArgs.Results.filenames;
            initStr = prsdArgs.Results.initStr;
            NumVec = prsdArgs.Results.NumVec;
            reProcessing = prsdArgs.Results.reProcessing;
            fixFalsePositive = prsdArgs.Results.fixFalsePositive;
            keepIlength = prsdArgs.Results.keepIlength;
            intpI = prsdArgs.Results.intpI;
            
            % work-out the filenames
            [filenames, data_path] = postprocess_dyfiles.get_files({'filenames',filenames,'initStr',initStr,'NumVec',NumVec});
            
            h = waitbar(0,'Please wait...');
            for j=1:length(filenames)                
                
                close all
                clear meas processed_meas
                
                filename=char(filenames(j));
                filenameWithPath=[data_path filename '.mat'];
                
                waitbar(j/length(filenames),h,['File: ' filename])

                tmp1 = who('meas','-file',filenameWithPath);
                tmp2 = who('processed_meas','-file',filenameWithPath);
                if isempty(tmp1)
                    disp(['File ' filename ' does not exists']); continue
                elseif ~isempty(tmp2) && ~reProcessing
                    disp(['File ' filename ' was already post processed']); continue
                end
                
                try
                    tmp = load(filenameWithPath,'meas');
                    meas = postprocess_dyfiles.fix_dy_file(tmp.meas);
                    meas = postprocess_dyfiles.remove_spikes(meas,fixFalsePositive,keepIlength);
                    
                    % perform Fourier transform to get the energy spectrum
                    [base_current,alpha1,real_sig,imag_sig,~,E0] = extract_pol_dyfiles(meas,intpI);
                    [~,energy,~,corrected_spectrum,~,~]=...
                        reconstruct_spectra(base_current,real_sig,imag_sig,E0,alpha1,base_current(end),0,1,0);
                    if energy(1)>energy(end)
                        meas.Energ_meV=flip(energy-E0); meas.SKw=flip(corrected_spectrum);
                    else
                        meas.Energ_meV=energy-E0; meas.SKw=corrected_spectrum;
                    end
                    
                    processed_meas = meas;
                    save(filenameWithPath,'processed_meas','-append');
                catch ME
                    disp(['File ' filename ' is nonexist/corrupt/lackOfInfo'])
                    continue
                end
                
            end
            close(h);
        end
        
        function meas = fix_dy_file(meas)
            
            % fix temperature reading
            tSample = floor(str2num(meas.loop(1).startStatus.tSample));
            if isempty(tSample) & length(meas.loop)>1
                tSample = floor(str2num(meas.loop(2).startStatus.tSample));
            end
            meas.temperature = tSample;
            
            if meas.numloops > length(meas.loop)% measurement ended before completion
                meas.numloops = length(meas.loop);
                meas.endBkg = meas.loop(1).bkg;
                meas.all = [];
                
                if length(meas.loop) > 1
                    loopTime = etime(meas.loop(2).startStatus.time,meas.loop(1).startStatus.time);
                    endTime = datevec(addtodate(datenum(meas.loop(1).startStatus.time),loopTime,'second'));
                    meas.endStatus = meas.loop(1).startStatus;
                    meas.endStatus.time = endTime;
                end
            end
            
            % fix gamma reading
            gamma = meas.loop(1).startStatus.gamma;
            if isempty(gamma), gamma = meas.endStatus.gamma;
            elseif isempty(meas.endStatus.gamma), meas.endStatus.gamma=gamma; end
            gammaSpecular = meas.loop(1).startStatus.gammaSpecular;
            if isempty(gammaSpecular), gammaSpecular = meas.loop(2).startStatus.gammaSpecular; end
            meas.dK = dK_for_gammaOfSpec_Ei(meas.beam.E0,gamma-gammaSpecular);            
        end
        
        function meas = remove_spikes(meas,fixFalsePositive,keepIlength)
            
            % Remove spikes by comparing loops
            % TODO: don't replace with the median, but just recalculate the
            % median values with out the outliers. Also, use the STD as a
            % weight for the fit procedure.
            if meas.numloops > 1
                Preal_mat = []; Pimag_mat = [];
                for i=1:meas.numloops
                    Preal_mat = [Preal_mat;meas.loop(i).Preal];
                    Pimag_mat = [Pimag_mat;meas.loop(i).Pimag];
                end
                med_Preal = repmat(median(Preal_mat),size(Preal_mat,1),1); spike_Preal = abs(Preal_mat-med_Preal)>0.02;
                med_Pimag = repmat(median(Pimag_mat),size(Pimag_mat,1),1); spike_Pimag = abs(Pimag_mat-med_Pimag)>0.02;

    %             figure; plot(Preal_mat','r-'); hold on; plot(Preal_mat(spike_Preal),'r+')
    %             hold on; plot(Pimag_mat','b-'); hold on; plot(Pimag_mat(spike_Pimag),'b+')

                Preal_mat(spike_Preal)=med_Preal(spike_Preal);
                Pimag_mat(spike_Pimag)=med_Pimag(spike_Pimag);

                for i=1:meas.numloops
                    meas.loop(i).Preal = Preal_mat(i,:);
                    meas.loop(i).Pimag = Pimag_mat(i,:);
                    meas.loop(i).Pmag = sqrt(meas.loop(i).Pimag.^2+meas.loop(i).Preal.^2);
                end
            end
            
            % Remove spikes by neighbours
            % The idea here is to find outliers in each loop and exclude
            % them from the data. Two itteration are used. First,
            % find potential outliers automatically using the hampel
            % function, then, confirm manually. Following that, reassign
            % and calculate relevant fields such as Preal, Pimage, and Pmag.
            % TODO: It is not necessary to exclude points, better just to
            % make them as 'not for use' using a dedicated field. Then the
            % fitting can ignore them using the 'Exclude' field.
            
            maxstd = 3;
            cutoff_setime = 5; %in [ps]
            numOfNeighbours = 5;
            
            if isfield(meas,'setime')
                indx = ~((abs(meas.setime)) < cutoff_setime);
            else
                cutoff_ibase = 0.1; %in ampere
                indx = ~((abs(meas.ibase)) < cutoff_ibase);
            end
            setime = []; ibase = [];
            Preal = []; Pimag = [];
            deltaPhase = []; deltaI0 = [];
            if keepIlength
                ibase_orig = [];
            end
            
            for i=1:meas.numloops
                
                loop = meas.loop(i);
                
                % Auto spike removal                                           
                
                [tmpPimag,tmp_excld_frm_imag,medianPimag,sigmaPimag] = hampel(loop.Pimag(indx),numOfNeighbours,maxstd);
                excld_frm_imag = zeros(size(indx)); excld_frm_imag(indx) = tmp_excld_frm_imag;
                [tmpPreal,tmp_excld_frm_real,medianPreal,sigmaPreal] = hampel(loop.Preal(indx),numOfNeighbours,maxstd);
                excld_frm_real = zeros(size(indx)); excld_frm_real(indx) = tmp_excld_frm_real;
                                
                % Now, finish with manual spike removal (to refine the
                % automated process. This implicitly assume that the
                % auto-removal will always find the right ones, and
                % sometimes will have 'false-positive' (but never
                % false-negative). This is of-course not 100% true
                
                % Assume all spikes found are false positive
                if fixFalsePositive
                    excld_frm_real = zeros(size(excld_frm_real));
                    excld_frm_imag = zeros(size(excld_frm_imag));
                end
                
                if find(excld_frm_real)>0
                    excld_frm_real = CustomFuncs.peakPointsFromFigure('dataX',meas.ibase,'dataY',meas.loop(i).Preal,'chosenPointsIndx',find(excld_frm_real),'titleStr',['Real, loop #' num2str(i)]);
                end
                if find(excld_frm_imag)>0
                    excld_frm_imag = CustomFuncs.peakPointsFromFigure('dataX',meas.ibase,'dataY',meas.loop(i).Pimag,'chosenPointsIndx',find(excld_frm_imag),'titleStr',['Imaginary, loop #' num2str(i)]);
                end

                Preal = [Preal loop.Preal(~(excld_frm_real | excld_frm_imag))];
                Pimag = [Pimag loop.Pimag(~(excld_frm_real | excld_frm_imag))];
                deltaPhase = [deltaPhase loop.deltaPhase(~(excld_frm_real | excld_frm_imag))];
                deltaI0 = [deltaI0 loop.deltaI0(~(excld_frm_real | excld_frm_imag))];
                if isfield(meas,'setime'), setime = [setime meas.setime(~(excld_frm_real | excld_frm_imag))]; end
                ibase = [ibase meas.ibase(~(excld_frm_real | excld_frm_imag))];
                if keepIlength, ibase_orig = [ibase_orig meas.ibase]; end
                
            end
            
            [meas.ibase, meas.mean.Preal,meas.mean.Pimag,meas.mean.deltaPhase,meas.mean.deltaI0] = CustomFuncs.merge_similar_points(ibase,Preal,Pimag,deltaPhase,deltaI0);
            
            % if we would like to keep the length of ibase, Preal, Pimag,
            % deltaPhase, and deltaI0 unchanged, we can interpolate them.
            if keepIlength
                ibase_orig=uniquetol(ibase_orig);
                meas.mean.Preal=interp1(meas.ibase,meas.mean.Preal,ibase_orig,'linear','extrap');
                meas.mean.Pimag=interp1(meas.ibase,meas.mean.Pimag,ibase_orig,'linear','extrap');
                meas.mean.deltaPhase=interp1(meas.ibase,meas.mean.deltaPhase,ibase_orig,'linear','extrap');
                meas.mean.deltaI0=interp1(meas.ibase,meas.mean.deltaI0,ibase_orig,'linear','extrap');
                meas.ibase=ibase_orig;
            end
            
            meas.mean.Pmag = sqrt(meas.mean.Preal.^2+meas.mean.Pimag.^2);
            if isfield(meas,'setime'), [meas.setime] = CustomFuncs.merge_similar_points(setime); end      
        end
        
        function [filenames, data_path] = get_files(varargin)
            % GET_FILES function will construct a file list, and work-out
            % what is the path for these files. Its a convenient way to
            % provide sequence of data-sets and bind them together.
            % Inputs: 
            %        filenames - list of specific filenames, if one whants
            %        to provide it directly.
            %        initStr - initial string that preceds all the files.
            %        NumVec - vector containing the number of file(s).
            % Outputs:
            %        filenames - {filename1, filename2, ...}
            %        data_path - the path to the files, assuming its in the
            %        folder defined by postprocess_dyfiles.projFold
            
            %% parse varargin
            prsdArgs = inputParser;   % Create instance of inputParser class.
            prsdArgs.addParameter('filenames', {}, @iscell);
            prsdArgs.addParameter('initStr', 'dy0', @ischar);
            prsdArgs.addParameter('NumVec', [], @isnumeric);
            prsdArgs.parse(varargin{:}{:});

            filenames = prsdArgs.Results.filenames;
            initStr = prsdArgs.Results.initStr;
            NumVec = prsdArgs.Results.NumVec;
            data_path='';
%             data_path=CustomFuncs.find_data_path(postprocess_dyfiles.projFold);
            
            %% SET the files to re-process
            for i=1:length(NumVec)
                filenames = {filenames{:}, [initStr num2str(NumVec(i))]};
            end
        end
        
        function meas_new = merge_dyfiles(varargin)
            
            [filenames, data_path] = postprocess_dyfiles.get_files(varargin);
            ibase = []; Preal = []; Pimag = []; deltaPhase = []; deltaI0 = [];
            
            for j=1:length(filenames)
                filename=[data_path char(filenames(j)) '.mat'];
                try
                    tmp=load(filename,'processed_meas');
                    if isempty(fieldnames(tmp))                        
                        continue
                    end
                    meas = tmp.processed_meas;

                    % Prepare mean real/imag, and mean temperature for
                    % files that where ended w/o completion of all loops.
                    for i=1:length(meas.loop), meanTmp(i,:)=meas.loop(i).Preal; end; Preal=[Preal mean(meanTmp,1)];
                    for i=1:length(meas.loop), meanTmp(i,:)=meas.loop(i).Pimag; end; Pimag=[Pimag mean(meanTmp,1)];
                    ibase = [ibase meas.ibase];
                    deltaPhase = [deltaPhase meas.loop(1).deltaPhase];
                    deltaI0 = [deltaI0 meas.loop(1).deltaI0];
                catch ME
                    disp(['File ' filename ' is nonexist/corrupt/lackOfInfo'])
                    continue
                end
            end
            
            % average repeating points
            [ibase,indx] = sort(ibase);            
            deltaPhase = deltaPhase(indx); deltaI0 = deltaI0(indx);
            Preal = Preal(indx); Pimag = Pimag(indx);
            tmp = diff(ibase); indx = find(abs(tmp)<1E-12);
            ibase(indx) = ibase(indx+1);
            [uniq_ibase,uniq_indx] = unique(ibase);
            
            uniq_indx = [uniq_indx' length(ibase)+1];
            for i=1:length(uniq_ibase)
                new_Pimag(i) = mean(Pimag(uniq_indx(i):uniq_indx(i+1)-1));
                new_Preal(i) = mean(Preal(uniq_indx(i):uniq_indx(i+1)-1));
                new_deltaI0(i) = mean(deltaI0(uniq_indx(i):uniq_indx(i+1)-1));
                new_deltaPhase(i) = mean(deltaPhase(uniq_indx(i):uniq_indx(i+1)-1));
            end
            
            % Force to be symmetric
            if 1 & find(uniq_ibase<0,1)
                zeroInd = find(uniq_ibase==0);
                uniq_ibase(1:zeroInd-1) = -1*fliplr(uniq_ibase(zeroInd+1:end));
                new_Pimag_left = mean([new_Pimag(1:zeroInd);fliplr(new_Pimag(zeroInd:end))],1); new_Pimag = [new_Pimag_left fliplr(new_Pimag_left(1:end-1))];
                new_Preal_left = mean([new_Preal(1:zeroInd);fliplr(new_Preal(zeroInd:end))],1); new_Preal = [new_Preal_left fliplr(new_Preal_left(1:end-1))];
                new_deltaI0_left = mean([new_deltaI0(1:zeroInd);fliplr(new_deltaI0(zeroInd:end))],1); new_deltaI0 = [new_deltaI0_left fliplr(new_deltaI0_left(1:end-1))];
                new_deltaPhase_left = mean([new_deltaPhase(1:zeroInd);fliplr(new_deltaPhase(zeroInd:end))],1); new_deltaPhase = [new_deltaPhase_left fliplr(new_deltaPhase_left(1:end-1))];
            end
            
            % interpulate to the finest grid
            if 1 & length(unique(diff(ibase))) > 1
                dI = min(diff(uniq_ibase));
                ibase_fine = uniq_ibase(1):dI:uniq_ibase(end);
                new_Pimag = interp1(uniq_ibase,new_Pimag,ibase_fine);
                new_Preal = interp1(uniq_ibase,new_Preal,ibase_fine);
                new_deltaI0 = interp1(uniq_ibase,new_deltaI0,ibase_fine);
                new_deltaPhase = interp1(uniq_ibase,new_deltaPhase,ibase_fine);
                uniq_ibase = ibase_fine;
            end
            
            % Padding
            if 1 & length(unique(diff(ibase))) > 1                
                dI = min(diff(uniq_ibase));
                ibase_fine = -10:dI:10;
                indx1 = ibase_fine<min(uniq_ibase)-1e-12; indx2 = ibase_fine>max(uniq_ibase)+1e-12;                
                new_Pimag = [repmat(new_Pimag(1),1,length(find(indx1))) new_Pimag repmat(new_Pimag(end),1,length(find(indx2)))];
                new_Preal = [repmat(new_Preal(1),1,length(find(indx1))) new_Preal repmat(new_Preal(end),1,length(find(indx2)))];
                new_deltaI0 = [repmat(new_deltaI0(1),1,length(find(indx1))) new_deltaI0 repmat(new_deltaI0(end),1,length(find(indx2)))];
                new_deltaPhase = [repmat(new_deltaPhase(1),1,length(find(indx1))) new_deltaPhase repmat(new_deltaPhase(end),1,length(find(indx2)))];               
                uniq_ibase = ibase_fine;
            end
            
            meas_new = meas;
            meas_new.filenames = filenames;
            meas_new.mean.Preal = new_Preal;
            meas_new.mean.Pimag = new_Pimag;
            meas_new.loop = [];
            meas_new.loop(1).Preal = new_Preal;
            meas_new.loop(1).Pimag = new_Pimag;
            meas_new.loop(1).Pmag = sqrt(new_Preal.^2+new_Pimag.^2);
            meas_new.loop(1).deltaPhase = new_deltaPhase;
            meas_new.loop(1).deltaI0 = new_deltaI0;
            meas_new.numloops=1;
            meas_new.ibase = uniq_ibase;
        end
        
        function sort_by_field(sortField, varargin)
            [filenames, data_path] = postprocess_dyfiles.get_files(varargin);
            for j=1:length(filenames)
                filename=char(filenames(j));
                filenameWithPath=[data_path filename '.mat'];

                try
                    tmp = load(filenameWithPath);
                    meas = tmp.processed_meas;
                    tmp = getfield(meas,sortField);
                    switch sortField
                        case 'ibase'
                            f(j) = max(tmp);
                        otherwise
                            f(j) = tmp;
                    end
                    
                catch ME
                    disp(['File ' filename ' is nonexist/corrupt/lackOfInfo'])
                    continue
                end
            end
            
            [f,indx] = sort(f);
            filenames = filenames(indx);
            [unique_f, indx] = unique(f);
            for i=1:length(indx)
                a=indx(i);
                if length(indx)==i
                    b=length(filenames);
                else b=indx(i+1)-1; end
                disp([sortField '=' num2str(unique_f(i)) '  :  ' filenames(a:b)])
            end            
        end
    end
end
        
