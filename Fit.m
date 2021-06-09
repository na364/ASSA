classdef Fit
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(Constant)
    end
    
    methods(Static)
        
        function [vBest,yBest,ci,rmse] = lorentzianFit(x,yOrig)
            
            options = fitoptions('Method','NonlinearLeastSquares','Lower',[1e10 10e3 -50e3 5e5],'Upper',[1e16 1e6 50e3 2e7],'StartPoint',[9e13 70e3 0 5e6]);
            
            ffun=fittype('(2*A/pi)*(w/(4*(x-x0)^2+w^2))+y0','options',options);
            [cfun, gof] = fit(x',yOrig',ffun);
            ci = confint(cfun,0.66);
            yBest=(2*cfun.A/pi)*(cfun.w./(4*(x-cfun.x0).^2+cfun.w^2))+cfun.y0;
            vBest = [cfun.A,cfun.w,cfun.x0,cfun.y0];
            
            if gof.rsquare < 0.85
                h = figure; plot(x,yOrig,x,yBest)
                choice = questdlg(['Fit: bad fit, rsquare=' num2str(gof.rsquare) ' ... drop the data?'], '', 'yes','no','choose intensity value', 'yes');
                close(h);
                switch choice
                    case 'yes', vBest=0; yBest=0; rmse=0;
                    case 'no',
                    case 'choose intensity value', yBestVal = inputdlg('choose intensity value: ','',1,{'-1'}); vBest = [str2num(yBestVal{:}) 0 0 0]; rmse = 0; return
                end                
            end
            rmse = gof.rmse;
        end
        
        function [vBest,yBest,ci,rsquare, cfun] = nExponentDecayFit(x,yOrig, varargin)
            
            %% parse varargin
            prsdArgs = inputParser;   % Create instance of inputParser class.
            prsdArgs.addParameter('Lower', [0 0 0], @isnumeric);%
            prsdArgs.addParameter('Upper', [inf inf inf], @isnumeric);%
            prsdArgs.addParameter('StartPoint', [0.1 0.08 0], @isnumeric);%
            prsdArgs.addParameter('numExp', 1, @isnumeric);%
            prsdArgs.addParameter('withBkgrnd', 1, @isnumeric);%
            
            prsdArgs.parse(varargin{:});
            
            Lower = prsdArgs.Results.Lower;
            Upper = prsdArgs.Results.Upper;
            StartPoint = prsdArgs.Results.StartPoint;
            numExp = prsdArgs.Results.numExp;
            withBkgrnd = prsdArgs.Results.withBkgrnd;
            
            if length(Lower) ~= length(Upper) || length(Upper) ~= length(StartPoint)
                error('Boundary or start conditions mismatch')
            end
            
            % if mismatch between number of requested exponenets, to
            % boundary/start conditions, numExp is followed, and default
            % boundary/start conditions are enforced.
            if length(StartPoint) ~= (numExp * 2 + withBkgrnd)
                %warning('Boundary/Start Conditions are incompatible with numExp - enforcing defaults')
                Lower = zeros(1, (numExp * 2 + withBkgrnd));
                Upper = [repmat([inf inf],1,numExp),inf];
                StartPoint = [repmat([0.1 0.08],1,numExp),0];
            end
            
            %%
            
            options = fitoptions('Method','NonlinearLeastSquares','Lower',Lower,'Upper',Upper,'StartPoint',StartPoint);
            funString = '';
            for i=1:numExp
                funString = [funString '+a' num2str(i) '*exp(-b' num2str(i) '*x)'];
            end
            if withBkgrnd
                funString = [funString '+e'];
            end
            ffun=fittype(funString,'options',options);
            if ~iscolumn(x),x=x';end
            if ~iscolumn(yOrig),yOrig=yOrig';end
            [cfun, gof] = fit(x,yOrig,ffun);
            ci = confint(cfun,0.66);
            yBest=cfun(x);
            vBest = [];
            for i=1:numExp
                vBest = [vBest, cfun.(['a' num2str(i)]),cfun.(['b' num2str(i)])];
            end
            
            if withBkgrnd
                vBest = [vBest, cfun.e];
            end
            
            if 0 && gof.rsquare < 0.90
                h = figure; plot(x,yOrig,x,yBest)
                choice = questdlg(['Fit: bad fit, rsquare=' num2str(gof.rsquare) ' ... drop the data?'], '', 'yes','no','choose intensity value', 'yes');
                close(h);
                switch choice
                    case 'yes', vBest=0; yBest=0; rmse=0;
                    case 'no',
                    case 'choose intensity value', yBestVal = inputdlg('choose intensity value: ','',1,{'-1'}); vBest = [str2num(yBestVal{:}) 0 0]; rmse = 0; return
                end                
            end
            rmse = gof.rmse;
            rsquare = gof.rsquare;
        end
        
        function [vBest,yBest,ci,rmse] = linearFit(x,yOrig)
            
            options = fitoptions('Method','Nonlinear','Lower',[-1 -1],'Upper',[1 1],'StartPoint',[0 0],'maxIter',2000);
            %options = fitoptions('Method','NonlinearLeastSquares');
            ffun=fittype('a*x+b','options',options);
            [cfun, gof] = fit(x',yOrig',ffun);
            ci = confint(cfun,0.66);
            yBest=cfun.a*x+cfun.b;
            vBest = [cfun.a,cfun.b];
            
            if 0 & gof.rsquare < 0.1
                h = figure; plot(x,yOrig,x,yBest)
                choice = questdlg(['Fit: bad fit, rsquare=' num2str(gof.rsquare) ' ... drop the data?'], '', 'yes','no','choose intensity value', 'yes');
                close(h);
                switch choice
                    case 'yes', vBest=0; yBest=0; rmse=0;
                    case 'no',
                    case 'choose intensity value', yBestVal = inputdlg('choose intensity value: ','',1,{'-1'}); vBest = [str2num(yBestVal{:}) 0 0]; rmse = 0; return
                end                
            end
            rmse = gof.rmse;
        end

        function [vBest,yBest,ci,rmse] = OneMinusExponentDecayFit(x,yOrig,PreExponentFactor)
            
            funcExpr = ['abs(a.*(1-' num2str(PreExponentFactor) '*exp(-x/b))+c)'];
            if PreExponentFactor == 1 % T1 from CPMG
                options = fitoptions('Method','NonlinearLeastSquares','Lower',[5e6 0.2 1e5],'Upper',[1e8 5 2e6],'StartPoint',[1e7 1 1e6],'MaxFunEvals',10000,'MaxIter',10000);
            elseif PreExponentFactor == 2 % T1 from Inversion Recovery
                options = fitoptions('Method','NonlinearLeastSquares','Lower',[1e7 0.9 3e6],'Upper',[5e9 4 1e7],'StartPoint',[1.6e8 1.16 2e6],'MaxFunEvals',10000,'MaxIter',10000);
            end
            ffun=fittype(funcExpr,'options',options);
            [cfun, gof] = fit(x',yOrig',ffun);
            ci = confint(cfun,0.66);
            yBest=cfun.a*(1-PreExponentFactor*exp(-x/cfun.b))+cfun.c;
            vBest = [cfun.a,cfun.b,cfun.c];
            
            if gof.rsquare < 0.85
                h = figure; plot(x,yOrig,'b+',x,yBest,'g+'); legend('yOrig','yBest');
                choice = questdlg(['Fit: bad fit, rsquare=' num2str(gof.rsquare) ' ... drop the data?'], '', 'yes','no','choose intensity value', 'yes');
                close(h);
                switch choice
                    case 'yes', vBest=0; yBest=0; rmse=0;
                    case 'no',
                    case 'choose intensity value', yBestVal = inputdlg('choose intensity value: ','',1,{'-1'}); vBest = [str2num(yBestVal{:}) 0 0]; rmse = 0; return
                end                
            end
            rmse = gof.rmse;
        end

        function [f gof] = gaussiansFit(numOfGaussians,ParamsConds, E_fit,int_fit)
            aStartPoint=[ ParamsConds.aStartPoint ParamsConds.bStartPoint ParamsConds.cStartPoint ParamsConds.dStartPoint];
            lower = [   ParamsConds.aLower      ParamsConds.bLower      ParamsConds.cLower      ParamsConds.dLower   ];
            upper = [   ParamsConds.aUpper      ParamsConds.bUpper      ParamsConds.cUpper      ParamsConds.dUpper  ];            
            
            opts = fitoptions('method','nonLinearLeastSquares','robust','off', 'startpoint',aStartPoint, 'lower', lower, 'upper',upper , ...
                'maxiter',800,'MaxFunEvals',800,'algorithm','trust-region');
            funStr = '';
            for i=1:numOfGaussians
                funStr = [funStr 'a' num2str(i) '*exp( -1*(x-b' num2str(i) ')^2/c' num2str(i) '^2 ) +'];
            end
            funStr = [funStr num2str(numOfGaussians) '*d'];
            ftype = fittype(funStr,'options',opts);
            [f gof] = fit(E_fit,int_fit,ftype);
        end

        function [f gof] = displacedGaussiansFit(numOfGaussians,ParamsConds, E_fit,int_fit, gaussians_for_elastic_and_background)
            aStartPoint=[ ParamsConds.aStartPoint ParamsConds.bStartPoint ParamsConds.cStartPoint ParamsConds.dStartPoint];
            lower = [   ParamsConds.aLower      ParamsConds.bLower      ParamsConds.cLower      ParamsConds.dLower   ];
            upper = [   ParamsConds.aUpper      ParamsConds.bUpper      ParamsConds.cUpper      ParamsConds.dUpper  ];            
            
            opts = fitoptions('method','nonLinearLeastSquares','robust','off', 'startpoint',aStartPoint, 'lower', lower, 'upper',upper , ...
                'maxiter',800,'MaxFunEvals',800,'algorithm','trust-region');
            funStr = '';            
            for i=1:numOfGaussians
                bString = '';
                if ~isempty(find(i==gaussians_for_elastic_and_background, 1))
                    bString=['-b' num2str(i)];
                else
                    for j=1:i                    
                        bString = [bString '-b' num2str(j)];
                    end
                        
                end
                funStr = [funStr 'a' num2str(i) '*exp( -1*(x' bString ')^2/c' num2str(i) '^2 ) +'];
            end
            funStr = [funStr num2str(numOfGaussians) '*d'];
            ftype = fittype(funStr,'options',opts);
            [f gof] = fit(E_fit,int_fit,ftype);
        end
     
    end
    
end

