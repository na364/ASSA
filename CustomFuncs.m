classdef CustomFuncs
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods(Static)
        
        function [xrot,yrot]=rot2D(x,y,rotAngle)
            Rmat = [cosd(rotAngle) sind(rotAngle); -sind(rotAngle) cosd(rotAngle)];
            x1 = reshape(x,1,[]);
            y1 = reshape(y,1,[]);
            vec = [x1;y1];
            vec1 = Rmat*vec;
            xrot=reshape(vec1(1,:),size(x));
            yrot=reshape(vec1(2,:),size(y));
        end
        
        function prcntMem = calc_required_memory(params)
            
            total_N_elements = 0;
            for i=1:length(params.prtcl)
                N_elements = params.prtcl(i).Nprtcl * length(params.prtcl(i).A) * ceil(params.stop_time/params.sample_time);
                total_N_elements = total_N_elements + N_elements;
            end
            
            if ispc
            else
                [a,b] = system('cat /proc/meminfo |grep MemTotal |awk ''{print $2}''');
                mem = str2num(b); % in kB    
            end
            
            prcntMem = (total_N_elements * 8 / 1000) / mem;
        
        end
        
        function parsave(varargin)
            fname = varargin{1};
            for i = 2:nargin
                savevar.res = varargin{i}; % other input arguments
            end
            
            if ~exist(fname,'file')
                save(fname,'-struct','savevar');
            else
                save(fname,'-struct','savevar','-append');
            end
                
        end
        
        function [uniqX, varargout] = merge_similar_points(vecX,varargin)
            % merge_similar_points will find similar values in vecX, within
            % tollerance of 1E-12 (can be changed, TODO). It will then
            % average the data points at these indcies for each vector in
            % varargin.
            varargout = varargin;
            [uniqX,indxOfUniqVal,indxOfSimVal]=uniquetol(vecX,1e-12);
            for i=1:length(varargout)
                vecY = varargout{i};                
                for j=1:length(uniqX)
                    ind = find(j == indxOfSimVal);
                    if length(ind) < 1, continue; end
                    vecY(ind) = mean(vecY(ind));
                end
                varargout{i} = vecY(indxOfUniqVal);
            end
            
%             varargout = varargin;
%             [uniqVal,zvl,indxOfuniqVal]=uniquetol(vecX,1e-12);
%             for i=1:length(uniqVal)
%                 ind = find(uniqVal(i) == indxOfuniqVal);
%                 if length(ind) < 2, continue; end
%                 for j=1:length(varargout)
%                     vecY = varargout{j};
%                     vecY(ind) = mean(vecY(ind));                    
%                     varargout{j} = vecY;
%                 end
%                 vecY = varargout{j};
%                 varargout{j} = vecY(unique(indxOfuniqVal));
% 
%             end
            
        end
        
        function home_path=find_home_path()
            [ret, name] = system('hostname');
            name = strtrim(name);
            if strcmp(name,'navidor') || strcmp(name,'na364-Latitude-5480')
                home_path='/home/na364/mnt/surfserver/';
            elseif strcmp(name,'DESKTOP-U20HL20')
                home_path='C:\Users\navidor\na364\';
            elseif strcmp(name,'surfcalc1')
                home_path='/home/na364/';
            elseif strcmp(name,'tomo')
                home_path='/home/na364/na364_server/';
            else
                home_path='/home/na364/';
            end
        end
        
        function scripts_path=find_scripts_path(projStr)
            home_path=CustomFuncs.find_home_path();
            scripts_path=[home_path 'scripts/' projStr '/'];            
        end
        
        function data_path=find_data_path(projStr, varargin)
            if strcmp(projStr(1),'/')
                data_path = projStr;
            elseif ~isempty(varargin)
                data_path = [varargin{1} '/'];
            else
                home_path=CustomFuncs.find_home_path();
                data_path=[home_path 'data/' projStr '/'];            
            end
        end
        
        function reminder=modAbs(x,y)
            reminder=mod(x,y);
            isBigger=reminder>y/2;
            mirrorReimder=y-reminder;
            reminder(find(isBigger))=mirrorReimder(find(isBigger));
        end
        
        function chosenPoints = peakPointsFromFigure(varargin)
            % peakPointsFromFigure will allow the user to choose points
            % using the mouse from a figure. The figure can either be supplied or created
            % 'on-the-fly' using the input data (if supplied).
            % Left click - Add the nearest point to the pointer position to the chosen list.
            % Right click - Remove the nearest point to the pointer position form the chosen list.
            % Middle click - return
            
            % parse varargin
            prsdArgs = inputParser;   % Create instance of inputParser class.
            prsdArgs.addParamValue('dataX', [], @isnumeric);
            prsdArgs.addParamValue('dataY', [], @isnumeric);
            prsdArgs.addParamValue('chosenPointsIndx', [], @isnumeric);
            prsdArgs.addParamValue('figHandle', [], @isnumeric);
            prsdArgs.addParamValue('closefig', 0, @isnumeric);
            prsdArgs.addParamValue('titleStr', 'Exclude (left click), re-include(right click), finish (middle button)', @ischar);
            prsdArgs.parse(varargin{:});
            figHandle = prsdArgs.Results.figHandle;
            closefig = prsdArgs.Results.closefig;
            chosenPointsIndx = prsdArgs.Results.chosenPointsIndx;
            dataX = prsdArgs.Results.dataX;
            dataY = prsdArgs.Results.dataY;
            titleStr = prsdArgs.Results.titleStr;
            
            % Open a figure if needed (using supplied data). In that case,
            % also close it on return
            if isempty(figHandle)
                closefig=1;
                figHandle = figure;
                plot(dataX,dataY,'o');
            else
                closefig=0;
                [data]=CustomFuncs.extract_figure_data(figHandle);
                dataX = data.x; dataY = data.y;
            end            
            if isrow(dataX), dataX = dataX'; end
            if isrow(dataY), dataY = dataY'; end
            
            figure(figHandle)
            XScale=diff(get(gca,'XLim'));
            YScale=diff(get(gca,'YLim'));
            title(titleStr)
            i=length(find(chosenPointsIndx));
            hold on;
            for j=1:i
                if j==0, continue; end
                h_tmp(j) = plot(dataX(chosenPointsIndx(j)),dataY(chosenPointsIndx(j)),'ro');
            end
            figure(figHandle)
            while 1                
                [x, y, button] = ginput(1);
                if button == 1
                    i=i+1;
                    r=sqrt(((dataX-x)./XScale).^2+((dataY-y)./YScale).^2);
                    [zvl chosenPointsIndx(i)] = min(r);
                    h_tmp(i) = plot(dataX(chosenPointsIndx(i)),dataY(chosenPointsIndx(i)),'ro');
                elseif button == 3
                    r=sqrt(((dataX(chosenPointsIndx)-x)./XScale).^2+((dataY(chosenPointsIndx)-y)./YScale).^2);
                    [zvl pointToExcludeIndx] = min(r);
                    delete(h_tmp(pointToExcludeIndx))
                    chosenPointsIndx = chosenPointsIndx(chosenPointsIndx(:)~=chosenPointsIndx(pointToExcludeIndx));
                    h_tmp = h_tmp(1:length(h_tmp) ~= pointToExcludeIndx);
                    i = i-1;
                else
                    %chosenPointsIndx = chosenPointsIndx(chosenPointsIndx(:)~=chosenPointsIndx(tmpchosenPointsIndx));
                    break
                end
            end
            
            chosenPoints = zeros(1,length(dataX));
            chosenPoints(chosenPointsIndx) = 1;            
            
            if closefig, close(figHandle); end
        end
        
        function includeIndx = manuallyCleanNoisyDataPoints(dataX,dataY)
            chosenPointsIndx = CustomFuncs.peakPointsFromFigure('dataX',dataX,'dataY',dataY);
            includeIndx = 1:length(dataX);
            includeIndx = includeIndx(~chosenPointsIndx);
        end
        
        function [data]=extract_figure_data(h)
            axesObjs = get(h, 'Children'); %axes handles
            if length(axesObjs) > 1
                disp(axesObjs)
                inpt =input('which AxesObject are are you interested in?');
                if inpt == 0, inpt=[1:length(axesObjs)]; end
            else
                inpt = 1;
            end                
                
            for i=1:length(inpt)
                dataObjs = get(axesObjs(inpt(i)), 'Children'); %handles to low-level graphics objects in axes
                objTypes = get(dataObjs, 'Type'); %type of low-level graphics object
                
                for j=1:size(objTypes,1)
                    if iscell(objTypes)
                        objTypes_ = objTypes{j,:};
                    else
                        objTypes_ = objTypes(j,:);
                    end
                    switch objTypes_
                        case 'surface'
                            %open(filename);
                            data(i).x=dataObjs.XData;
                            data(i).y=dataObjs.YData;
                            data(i).z=dataObjs.ZData;
                            data(i).c=dataObjs.CData;
                        case 'line'
                            data(i).x = get(dataObjs, 'XData'); %data from low-level grahics objects
                            data(i).y = get(dataObjs, 'YData');
                            data(i).z = get(dataObjs, 'ZData');
                        case'hggroup'
                            if strcmp(objTypes{2},'hggroup')
                                for j=1:length(dataObjs)
                                    data.x = get(dataObjs, 'XData'); %data from low-level grahics objects
                                    data.y = get(dataObjs, 'YData');
                                    data.L = get(dataObjs, 'LData');
                                    data.U = get(dataObjs, 'UData');
                                end
                            end
                        case 'errorbar'
                            data(j).x = get(dataObjs(j), 'XData'); %data from low-level grahics objects
                            data(j).y = get(dataObjs(j), 'YData');
                        otherwise
                            %warning(['The objType ' objTypes ' is not defined.'])
                    end
                end
                
            end
        end
        
        function [n,m] = locatePointsInData(X,Y,x0,y0)
            
            bracket = abs(X(2)-X(1));
            n = zeros(size(x0)); m = zeros(size(x0));
            
            for i=1:length(x0)
                [nTmp,mTmp] = find(abs(X-x0(i)) < bracket/5 & abs(Y-y0(i)) < bracket/5,1);
                if isempty(nTmp)
                    %disp(['Zero points at i=' num2str(i) ' ... x0 = ' num2str(x0(i)) ' ... y0 = ' num2str(y0(i))]);
                    nTmp = 0; mTmp = 0;
                end
                n(i) = nTmp; m(i) = mTmp;
            end
        end
        
        function paramIndMat = allPerm(varargin)
            % ALLPERM spans space of variables
            % Recurcively span all the permutations of values in cell-array
            % arguments. For example, if we have 3 parameters with 3
            % variable each, var1={val1.1,val1.2,val1.3},
            % var2={val2.1,val2.2,val2.3}, etc, paramIndMat =
            % allPerm(var1,var2,var3) will generate paramIndMat were
            % {:,j} is a combination of var1-3. So size(paramIndMat)=3x27
            % Each cell-array of values can be of strings, numbers, but not
            % both (as its meaningless). For example, varargin could be
            % {{1,2},{'aa','bb'}} but not {{1,'bb'},{'aa',2}}

            if nargin == 1
                paramIndMat = varargin{1};
                return
            end

            % After all calls are done, paramIndMat holds (in every 'call')
            % the previous permutated parameters.
            paramIndMat = CustomFuncs.allPerm(varargin{1:end-1});
            L = size(paramIndMat,2);
            
            % concatenate repetitions of the previous permutated
            % parameters, as preparation to add the 'new' parameters
            paramIndMat = repmat(paramIndMat,1,length(varargin{nargin}));
            
            % repeat the new parameters in the 'amount' needed for all the
            % previous permutations
            tmp_new_params = repmat(varargin{nargin},1,L);
            if iscell(tmp_new_params) && isnumeric(tmp_new_params{1})
                tmp_new_params = cell2mat(tmp_new_params);
            end
            
            tmp = sort(tmp_new_params);
            if isnumeric(tmp), tmp = num2cell(tmp); end
            if isnumeric(paramIndMat), paramIndMat = num2cell(paramIndMat); end
            paramIndMat(nargin,:) = tmp;
        end
    
        function a = mvdlg(h)
            % display on figure(h) an OK and Cancel buttons, and return the
            % input from the user
            dialogWind = h;
            okayButton = uicontrol('style','pushbutton','Units','normalized',...
                'position', [.05,.9,.1,.075],'string','OK','callback',@okCallback);
            cancelButton = uicontrol('style','pushbutton','Units','normalized',...
                'position', [.2,.9,.1,.075],'string','Cancel','callback',@cancCallback);
            a = {};
            uiwait(dialogWind);
            function okCallback(hObject,eventdata)
                a = 'Yes';
                uiresume(dialogWind);
            end
            
            function cancCallback(hObject,eventdata)
                a = 'No';
                uiresume(dialogWind);
            end
        end

    end
    
end

