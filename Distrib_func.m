classdef Distrib_func
    % DISTRIB_FUNC is the library of functions which calculate a
    % distribution function of currents and generate the ibase matrix
    % accordingly
    
     properties(Constant)
     end
     
     methods(Static) 
         
            function [nfinal, pfinal] = gen_rand_current(varargin)
                %   generates a vector of np unique integers selected
                %   according to the discrete distribution function pi.  
                %   The procedure makes use of the rejection method 
                %   (Numerical recipes Sect7.3.6) in its simple form.
                %   Outputs:
                %       nfinal: the list of random integers from the distribution       
                %       pf: the distribution function with pf(n)=0 for all the selected
                %           random integers (provided so that further calls can be made
                %           to give additional unique integers from the same initial
                %           distribution).
                %       ibasetoreconstruct: the basis of currents to
                %       reconstruct, in order to ensure that I = 0 A is
                %       taken
                %
                % First check that there are sufficient numbers left (i.e. with pi>0) for
                % the routine to work.  (NB the routine will be very inefficient if there
                % are a large number of zero values in pinit)
                for i=1:2:length(varargin)
                    switch varargin{i}  
                        case 'pinit'
                            pinit = varargin{i+1};
                        case 'np'
                            np = varargin{i+1};
                        case 'ibasetoreconstruct'
                            ibasetoreconstruct = varargin{i+1};
                        otherwise
                            warning(['variable ' varargin{i} ' is not defined']);
                    end
                end
                
                if exist('ibasetoreconstruct','var')==1 
                    len=length(pinit);
                    if (len>length(ibasetoreconstruct))
                        len = length(ibasetoreconstruct);
                    else
                        len=length(pinit);
                    end
                else
                      len=length(pinit);
                end
                nd=sum(abs(pinit) > eps);
                if np > nd
                    error('length of sequence requested is longer than the distribution')
                end
                nfinal=zeros(np,1);
                pfinal=pinit;
                pmax=max(pfinal);
                %
                %   loop through to create the np random integers
                %
                ncount=1;
                while ncount <= np
                   %generate 2 random numbers for y and x axes
                   ry = rand()*pmax;   % random no. between 0 and pmax
                   rx = fix(rand()*len)+1;
                   if pfinal(rx) > ry
                       pfinal(rx)=0.0;
                       pmax=max(pfinal);
                       nfinal(ncount)=rx;
                       ncount = ncount + 1;
                   end
                end
                
                % check that I = 0 A is sampled:
                if exist('ibasetoreconstruct','var')==1
                    ibasetomeasure = ibasetoreconstruct(nfinal);
                    if isempty(find(ibasetomeasure == 0.0))
                        nfinal(end+1) = find( ibasetoreconstruct==0);
                    end
                end
            end

            function [pfact, pdf_continuous,pdf_discrete] = pdf_GUI(varargin)
                %This function produces both the continuous and discrete lorentzian pdfs for 
                %the user interface. Only the positive current region is represented in the 
                %distribution for basic user guidance. This means the pdf is generated and
                %plotted for the region between 0 and maximum current input via the use of 
                %pfact array. 

                %Conversion from current I, to index, n using: 
                %n=[((I-Imin)/(Imax-Imin))*(N-1)]+1
                
                % Imin=1;
                % Imax=8;
                % dI=0.0001;
                % r=[1 1 2 3 3 3 4 4 4 4];
                  for i=1:2:length(varargin)
                    switch varargin{i}  
                        case 'Imin'
                            Imin = varargin{i+1};
                        case 'Imax'
                            Imax = varargin{i+1};
                        case 'dI'
                            dI = varargin{i+1};                         
                        case 'DistWidth'
                            DistWidth = varargin{i+1}; 
                        case 'r'
                            r = varargin{i+1};                         
                        otherwise
                            warning(['variable ' varargin{i} ' is not defined']);
                    end
                end

                % correct for Imin>Imax  so that the scan direction is always positive    

                if Imin>Imax 
                    Itemp=Imin;
                    Imin=Imax;
                    Imax=Itemp;
                end         

                % no of buttons
                nb=length(r);

                %This ensures : O-->  n=1 && Imin -->  n=min  && Imax -->  n=nmax

                imin=min(abs(Imin),abs(Imax));
                imax=max(abs(Imin),abs(Imax));

                % Need to make sure the length of index array is divisible by the number of
                % regions (i.e. number of buttons) 

                nmin=(ceil(imin/(dI*nb))*nb);
                nmax=(ceil(imax/(dI*nb))*nb);
                %DistWidth=imax/2;

                %pfact array initially set to include all current values 
                pfact=ones(1,nmax);


                if (Imin*Imax)>=0 

                %set pfact=0 for all values outside the requested range (i.e. below imin)
                %as this condition implies single polarity meaning zero current needs to be
                %reinforced into the scan range

                      pfact(1:(nmin-1))=0;
                      pf=find(pfact~=0);

                      nw=(DistWidth/imax)*length(pf);
                      n_pdf=1:nmax;
                      pdf_continuous= ((nw^2)./(nw^2 +n_pdf.^2));

                elseif (Imin*Imax)<0 

                    nw=(DistWidth/imax)*nmax;
                    n_pdf=1:nmax;
                    pdf_continuous= ((nw^2)./(nw^2 +n_pdf.^2));
                end

                %The following lines of code concerns with the generation of discrete
                %lorentzian. This is achieved by equally dividing the current grid into 10
                %sections and averaging over the specific regions as determined by the
                %region vector (button set-up).

                y=(reshape(pdf_continuous,[],nb))';

                [a] = histc(r,unique(r));

                pdf_discrete=[];

                for j=1:numel(a)

                        i=a(j);
                        if i>1
                        p=y(1:i,:);
                        Np=numel(p);
                        Pd=sum(sum(p(1:i,:)))/Np;
                        pdf_discrete=[pdf_discrete,(zeros(1,Np)+Pd)];
                        y(1:i,:)=[];
                        elseif i==1
                        p=y(i,:);
                        Np=numel(p);
                        Pd=sum(p(i,:))/Np;
                        pdf_discrete=[pdf_discrete,(zeros(1,Np)+Pd)]; %#ok<*AGROW>
                        y(i,:)=[];
                        end

                end


                % 
                % figure(1)
                % plot(pdf_continuous.*pfact,'b--o');
                % hold on;
                % stairs(pdf_discrete.*pfact,'r--o');


            end

            function [pdflist_continuous,pdflist_discrete,no] = pdf_list(varargin)

                    %This function generates 2 lorentzian pdfs,continuous and discrete to be 
                    %passed as an index list but corrects for the positive current
                    %reinforcement. This is achieved by means of flipping so that the current
                    %ranges that were initially excluded from the distribution can now be covered.
                    %The function also notes down the position of zero current,no, within the 
                    %index list.


                    % Imin=2;
                    % Imax=8;
                    % dI=0.0001;  
                    % r=[1 1 2 3 3 3 4 4 4 4];
                    
                   for i=1:2:length(varargin)
                    switch varargin{i}  
                        case 'Imin'
                            Imin = varargin{i+1};
                        case 'Imax'
                            Imax = varargin{i+1};
                        case 'dI'
                            dI = varargin{i+1};                         
                        case 'DistWidth'
                            DistWidth = varargin{i+1}; 
                        case 'r'
                            r = varargin{i+1};                          
                        otherwise
                            warning(['variable ' varargin{i} ' is not defined']);
                    end
                end
                 



                    % correct for Imin>Imax  as before to ensure positive slope    

                    if Imin>Imax 
                        Itemp=Imin;
                        Imin=Imax;
                        Imax=Itemp;
                    end         

                    % no of buttons
                    nb=length(r);



                    imin=min(abs(Imin),abs(Imax));
                    imax=max(abs(Imin),abs(Imax));

                    %Now the index n serves as a flag for flipping
                    nmin=(ceil(imin/(dI*nb))*nb);
                    nmax=(ceil(imax/(dI*nb))*nb);


                    %DistWidth=(imax)/2;

                    %Calls the pdf_GUI function which generates the initial lorentzian pdfs
                    [pfact, pdf_cont,pdf_disc] = Distrib_func.pdf_GUI('Imin',Imin,'Imax',Imax,'dI',dI,'DistWidth',DistWidth,'r',r);
                    
                    %for (Imin*Imax)>=0 where the pfact array was needed

                    if (Imin*Imax)>=0 && abs(Imax)>abs(Imin)

                           pdf_continuous=pdf_cont.*pfact;
                           pdf_discrete=pdf_disc.*pfact;
                           pdflist_continuous=(pdf_continuous(1:nmax));
                           pdflist_discrete=(pdf_discrete(1:nmax));
                           no=1;

                    elseif(Imin*Imax)>=0 && abs(Imax)<abs(Imin)

                           pdf_continuous=pdf_cont.*pfact;
                           pdf_discrete=pdf_disc.*pfact;
                           pdflist_continuous=(fliplr(pdf_continuous(1:nmax)));
                           pdflist_discrete=(fliplr(pdf_discrete(1:nmax)));
                           no=nmax; 
                    end

                    %for (Imin*Imax)<0 

                     pdf_continuous=pdf_cont;
                     pdf_discrete=pdf_disc;

                    if (Imin*Imax)<0 && abs(Imax)>abs(Imin)   

                        pdflist_continuous=[fliplr(pdf_continuous(2:nmin)) pdf_continuous(1:nmax)];
                        pdflist_discrete=[fliplr(pdf_discrete(2:nmin)) pdf_discrete(1:nmax)];
                        no=nmin;

                    elseif (Imin*Imax)<0 && abs(Imax)<abs(Imin)    

                        pdflist_continuous=[fliplr(pdf_continuous(2:nmax)) pdf_continuous(1:nmin)];
                        pdflist_discrete=[fliplr(pdf_discrete(2:nmax)) pdf_discrete(1:nmin)];
                        no=nmax;

                    elseif (Imin*Imax)<0 && abs(Imax)==abs(Imin)    

                        pdflist_continuous=[fliplr(pdf_continuous(2:nmax)) pdf_continuous];
                        pdflist_discrete=[fliplr(pdf_discrete(2:nmax)) pdf_discrete];
                        no=nmax;   

                    end

                    % 
                    % figure(1)
                    % plot(pdflist_continuous);
                    % hold on;
                    % stairs(pdflist_discrete,'r--o');

            end
                                   
            function [bincounts,binranges_current]=ListHistogram(list,levels,ibase,figurenum,sub,Shift,Title)
                % Contructs the corresponding histogram for a sampling map
                % generated using the subsampleinlevels function
                % list = sampling map
                % levels = sampling levels

                res=length(list);
                L=length(levels);
                if nargin<5
                    Shift=0;
                end
                binranges=[10^(-8),levels+10^(-8)];
                binranges_current=min(ibase)+binranges*(max(ibase)-min(ibase))/max(binranges);
                binwidths=ceil(binranges(2:L+1))-ceil(binranges(1:L));
                if size(list,1) ==1
                    list =list.';
                end    
                listnum=(1:res).*(list.');
                listnum(list==0)=[];
                bincounts=histc(listnum,binranges);
                bincounts=bincounts(1:L);
                bincounts=bincounts./binwidths;
                bincounts(L+1)=0;
                if (nargin>3)
                    figure(figurenum)
                end
                if(nargin>4)
                    subplot(sub(1),sub(2),sub(3));
                end
                if nargin>5
                    binranges_current=binranges_current-Shift;
                    bar(binranges_current,bincounts,'histc')
                    %axis([1-Shift,res-Shift,0,1]);
                    axis([min(binranges_current),max(binranges_current),0,1])

                else
                    bar(binranges_current,bincounts,'histc')
                    axis([min(binranges_current),max(binranges_current),0,1])
                end
                if nargin>6
                    title(Title)
                else
                    title('Sampling Histogram')
                end
                %output=bincounts;
            end
                   
            function output=subsampleinlevels(res,levels,subper)
                 % Used to create the random samples
                    % res = base resolution level=maximum number of samples
                    % levels = a vector of numbers ranging between 1 and res used to break up 
                    %          the sampling domain into regions
                    % subper = a vector of fractions denoting the number of samples taken in 
                    %          each level 
                    % The user must now keep track of what the indexes correspond to. For
                    % example if we running over frequencies from kappa_min to kappa_max then
                    % output(i) coresponds to kappa_min +(i-1)*deltakappa.

                    run=1:res;
                    values=abs(run);

                    samplingtotal=0;
                    output=zeros(res,1);
                    for r=1:length(levels)
                        levelcheck=zeros(res,1);
                        if r~=1
                            levelcheck( (values <= levels(r)) & (values > levels(r-1)))=1;
                        else
                            levelcheck(values <= levels(r))=1;
                        end
                        range=sum(levelcheck);
                        sampling=floor(range*subper(r));
                        randindex=randperm(range);
                        randindex=randindex(1:sampling);
                        trueindex=zeros(range,1);
                        trueindex(randindex)=1;
                        if r~=1
                            output((values <= levels(r)) & (values > levels(r-1)))=trueindex;
                        else
                            output(values <= levels(r))=trueindex;
                        end
                        samplingtotal=samplingtotal+sampling;
                    end
                    
                    display(samplingtotal)
            end
           
            function [ibase_matrix,CS_s] = gen_ibase_matrix(varargin)
                % GEN_IBASE_MATRIX is a function which generates the
                % current vector to measure in the subsequent loops and the
                % structure with the info needed for the Compressed sensing
                % reconstruction.
                %
                %   Inputs:
                %       Imin:                      minimum current in (A)
                %       Imax:                      maximum current in (A)
                %       deltaI:                    interspace between
                %                                  currents in (A)
                %       numloops:                  number of loops
                %       histogramstr               selection of currents by
                %                                  the means of a
                %                                  HISTOGRAM. choose 'on'
                %                                  or 'off' 
                %         =>  levels:                    bins of currents
                %         =>  prob_levels:               probability in each bin
                %       pdf_constr                 selection of currents by
                %                                  means of a CONTINUOUS lorentzian distribution
                %                                  function. choose 'on' or
                %                                  'off'
                %       => FWYHM                   FWHM of the distribution
                %       => num_levels              bins to construct an
                %                                  equivalent histogram
                %       pdf_discrtstr              selection of currents by
                %                                  means of a DISCRTETE lorentzian distribution
                %                                  function. choose 'on' or
                %                                  'off'
                %       => FWYHM                   FWHM of the distribution
                %       => num_levels              bins to construct an
                %                                  equivalent histogram
                %
                %      => npoints_loop             number of currents to
                %                                  scan per loop
                %      => fractp_loop             alternatively one can
                %                                 specify  the fraction between 0 and 1 of points per loop and the
                %                                  total npoints_loop = fractp_loop*resolution
                %                                  scan per loop
                %
                %       Method:                    method for the CS reconstruction. Choose among 'DFT', 'DWT' (discrete wavelets) and 'Continuous'
                %
                %       file_test:                 in the case we want to
                %                                  run the experiment in test mode with an already
                %                                   existing measurement
                %
                % The outputs are:
                %      ibase_matrix                matrix containing the
                %                                  different ibase to measure in each loop
                %      CS_s                        the structure containing
                %                                  the information to do compressed sensing
                %
                %       CS_s.ibasetoreconstruct    the ibase vector to reconstruct
                %
                %       CS_s.ibase_matrix          the ibase matrix to measure (the number of columns corresonds to the number of  loops and the number of rows is the number of
                %                                  current per loop)
                %
                %       CS_s.list_mask_matrix      the "mask" matrix: the number of columns corresponds to the number of
                %                                  loops and the number of rows corresponds to the
                %                                  number of elements of ibasetoreconstruct. It is
                %                                  made of 1s and 0s which label the measured/non
                %                                  measured current for each loop
                %
                %       CS_s.list_index_matrix     matrix which contains the index of the currents measured in each loop.
                %                                  Number columns = numloops; Number rows: number of currents per loop
                %      CS_s.Method                 method for CS reconstruction. 'DFT', 'DWT' or 'Continuous'
                
                % READ IMPUT PARAMETERS:
                 for i=1:2:length(varargin)
                    switch varargin{i}  
                        case 'Imin'
                            Imin = varargin{i+1};
                        case 'Imax'
                            Imax = varargin{i+1};
                        case 'deltaI'
                            deltaI = varargin{i+1};
                        case 'numloops'
                            numloops = varargin{i+1};   
                        case 'histogramstr'
                            histogramstr = varargin{i+1};   
                        case 'levels'
                            levels = varargin{i+1};
                        case 'prob_levels'
                            prob_levels = varargin{i+1}; 
                        case 'pdf_contstr'
                            pdf_contstr = varargin{i+1};  
                        case 'pdf_discrstr'
                            pdf_discrstr = varargin{i+1};  
                        case 'FWHM'
                            FWHM = varargin{i+1}; 
                        case 'npoints_loop'
                            npoints_loop = varargin{i+1}; 
                        case 'fractp_loop'
                            fractp_loop = varargin{i+1}; 
                        case 'num_levels'
                            num_levels = varargin{i+1}; 
                        case 'Method'
                            Method = varargin{i+1};  
                        case 'file_test'
                            file_test = varargin{i+1};  
                        otherwise
                            warning(['variable ' varargin{i} ' is not defined']);
                    end
                 end
                
                 
                 % SET DEFAULT VALUES:
                 
                  
                 
                   if exist('Imin','var')==0
                        Imin = -1.0;
                        display('Set default Imin = -1.0 A');
                   end
                   
                   if exist('Imax','var')==0
                        Imax = 1.0;
                        display('Set default Imax = -1.0 A');
                   end
                    
                   if exist('deltaI','var')==0                     
                        deltaI = 1e-3;
                        display(['Set default deltaI = 1e-3 A']);
                   end
                
                    
                   if exist('numloops','var')==0
                        numloops = 1;
                        display('Set default numloops = 1');
                   end
                   
                    if exist('histogramstr','var')==0 && exist('pdf_contstr','var')==0 && exist('pdf_discrstr','var')==0
                       prob_levels = 1;
                       levels = 1;
                       histogramstr = 'on';                       
                       display('Set default Full sampling with histogram method');
                    end
                                      
                   if exist('levels','var')==0
                        levels = [1/5,2/5,3/5,4/5,5/5];
                        levels_current = Imin + (Imax-Imin).*levels;
                        display(['Set default levels of current = ' levels_current 'A']);
                   end
                    
                   if exist('prob_levels','var')==0 
                        prob_levels = [1,1,1,1,1];
                        display(['Set default probability of current = ' prob_levels 'A']);
                   end
                   
                   if exist('FWHM','var')==0
                        FWHM = 1;
                        %display(['Set default probability of current = ' prob_levels 'A']);
                   end
                   
                    if exist('num_levels','var')==0
                        num_levels = 1;
                        %display(['Set default probability of current = ' prob_levels 'A']);
                   end
                   
                   if exist('Method','var')==0
                       Method = 'DFT';
                       display('Set default method for CS reconstruction: DFT ');
                   end
                   
                    if exist('file_test','var')==1
                        test_mode = 1;
                        meas_test1 = load(file_test);                   
                        meas = meas_test1.meas;
                        
                        Imin =min(meas.ibase);
                        Imax =max(meas.ibase);
                        display(['Test mode is on and will apply to file ' file_test]);
                    else 
                        test_mode = 0;
                    end
                   
                      % Calculate ibase to reconstruct or read it from file
                   if (test_mode == 0)
                     ibase_reconstruct = [ Imin:deltaI:-deltaI 0:deltaI:Imax];
                   else
                       ibase_reconstruct = meas.ibase;
                   end
                   res = length(ibase_reconstruct);
                   
                    if exist('pdf_contstr','var')==1 || exist('pdf_discrstr','var')==1
                         if exist('npoints_loop','var')==0 && exist('fractp_loop','var')==0
                             npoints_loop = res;
                             disp('Number of currents per loop is set by default to the maximum');
                         elseif exist('npoints_loop','var')==0 && exist('fractp_loop','var')==1
                             npoints_loop = round(res*fractp_loop);
                             disp(['Number of currents per loop is set to ' num2str(npoints_loop)]);
                         end                                                                   
                   end 
                  

                   
                   % build the histogram x-axis by converting the
                   % num_levels vector (as provided in the interface) into
                   % a binrange vector.
                   if exist('num_levels','var')==1 && exist('histogramstr','var')==0
                         index_uniq = unique(num_levels);
                         redundant_index = histc(num_levels,index_uniq);
                         levels = zeros(1,length(index_uniq));
                       
                         levels(1) = redundant_index(1)/length(index_uniq);
                         
                        for jj = 2:length(index_uniq)
                            levels(jj) = levels(jj-1) + redundant_index(jj)/length(index_uniq); 
                        end
                   end
                   
                 
                   
                   levels_current = Imin + (Imax-Imin).*levels;
                   
                   % inizialize the list of 0 and 1 which determine the
                   % measured currents:
                   list = zeros(res,1);
                   
                   % figure to plot the histogram of sample current for
                   % each loop:
                   figure(10);
                   
                   % for each loop sub sample ibase_reconstruct to obtain
                   % the ibase_matrix to measure:
                   
                   for kk = 1:numloops
                       
                       % genereates the list of currents to measure with a
                       % histogram, or a pdf (continuous/discrete)
                       if exist('histogramstr','var')==1
                            list= Distrib_func.subsampleinlevels(res,res*levels,prob_levels);
                            nfinal = find(list == 1);
                       elseif exist('pdf_contstr','var')==1        
                           if (test_mode == 1)
                                [pdf_cont,~,~] = Distrib_func.pdf_list('Imin',Imin,'Imax',Imax,'dI',deltaI,'DistWidth',FWHM,'r',num_levels,'file_test',file_test);
                           else
                               [pdf_cont,~,~] = Distrib_func.pdf_list('Imin',Imin,'Imax',Imax,'dI',deltaI,'DistWidth',FWHM,'r',num_levels);
                           end
                           [nfinal,~] = Distrib_func.gen_rand_current('pinit',pdf_cont,'np',npoints_loop,'ibasetoreconstruct',ibase_reconstruct);                              
                           list( nfinal) = 1; 
                       elseif exist('pdf_discrstr','var')==1
                            if (test_mode == 1)
                               [~,pdf_discrt,~] = Distrib_func.pdf_list('Imin',Imin,'Imax',Imax,'dI',deltaI,'DistWidth',FWHM,'r',num_levels,'file_test',file_test);
                           else
                               [~,pdf_discrt,~] = Distrib_func.pdf_list('Imin',Imin,'Imax',Imax,'dI',deltaI,'DistWidth',FWHM,'r',num_levels);
                           end
                           
                           [nfinal,~] = Distrib_func.gen_rand_current('pinit',pdf_discrt,'np',npoints_loop,'ibasetoreconstruct',ibase_reconstruct);
                           list( nfinal) = 1; 
                       end
                       
                       % stores the list of 1 (measured currents) and 0  for each loop
                       % within a matrix, to allow Compressed sensing
                       % reconstruction
                       list_mask_matrix(:,kk) = list;
                       
                        % stores the list of current  for each loop
                       % within a matrix, to allow the recombination of the
                       % measured data in each loop
                       list_index_matrix(:,kk) = nfinal;
                       
                       % stores in a matrix the currents to measure for each
                       % loop
                       ibase_matrix(:,kk) = ibase_reconstruct(nfinal);
                       
                       % Calculates the shift to represent the current
                       % histogram with the bars centered in the bin
                       if length(levels_current) > 1
                            Shift = (levels_current(2)-levels_current(1))/2;
                       else
                           Shift = 0;
                       end    
                        Shift = 0;
                      
                       
                       % Plots the histogram of the measured currents:
                       Distrib_func.ListHistogram(list,res*levels,ibase_reconstruct,10,[ceil(sqrt(numloops)) ceil(sqrt(numloops)) kk],Shift,['Histogram for a subsampling of ' num2str(sum(list)) ' / ' num2str(res)]) ;
                   end

                   CS_s.ibasetoreconstruct = ibase_reconstruct;
                   CS_s.ibase_matrix = ibase_matrix;
                   CS_s.list_mask_matrix = list_mask_matrix;
                   CS_s.list_index_matrix = list_index_matrix;
                   CS_s.Method = Method;
            end
            
     end
end