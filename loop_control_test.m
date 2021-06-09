function [config_struct,flag_end_scan,analysis]=loop_control(config_struct,nlist,loop_now,figure_struct,CS_s)%,meas)
% SIMPLE_GUI_SE_meas  performs what we call an
% experiment: a collection of scans of current at a given temperature,
% coverage and scattering angle which is stored in a single data file.
% It consists in: 
% 1) a preamble where the graphical interphase is designed
% 2) an external loop which
% is the experiment, 
% 3) within the loop we call "scan_current" to measure
% one single scan, 
% 4) we save each scan within a file, 
% 5) each Npoints we  calculate the "scattering function" and the "intermediate scattering function"
% by the means of the Compress Sensing (CS) method developed by A. Jones.

 
 
 %addpath /home/shared/SpinEcho/Spinecho2Icon/data/Ni111_Graphene_H2O/
    
    % example of phonon measurements, i.e. oscillations
    filelist={'dy011880.mat'}; 
    
    % with this single decay it fails!!! 
    %filelist={'dy010562.mat'}; 
     %filelist={'dy012858.mat'}; 
    
    % example of a single decay
    %filelist={'dy014926.mat'};
    
    
    
    
    filein = char(filelist);  
    
    
    %Read in the file (as a meas structure)
    if isnumeric(filein)
       fileindex=char(['dy0',num2str(filein,'%05d')]);
       filename=[fileindex '.mat'];
    else
        if ~strfind(filein,'dy')
            error('File not recognised, enter in numeric form i.e 1449 or dy file path as a string')
        else
            if ~strfind(filein,'.mat')
                filename=[filein '.mat'];
            else
                filename=filein;
            end
        end
    end
    load(filename);
    meas_1 = meas;
         


 
   %% =================== Calculate the current and prepare the list and the structure for the Compress sensing: ===================
   
    [ibase,list_meas,ibase_reconstruct] = index2current(config_struct);
   
    
   field1 = 'list'; value1 = list_meas;
   field2 = 'res'; value2 = length(ibase_reconstruct);
   field3 = 'Npoints'; value3 = 50; % the number of experimental points after which we want to "update" the SF and the ISF
   field4 = 'save'; value4 = 0; % flag =1 means create a file with the CS result, flag=0 do not create the file
   field5 = 'nlist'; value5 = nlist; % stores the list of current indexes which are measured
   field6 = 'ibasetoreconstruct'; value6 = ibase_reconstruct; % current vector to reconstruct
   field7 = 'Method'; value7 = config_struct.Method; % current vector to reconstruct
   field8 = 'InterpolationMultiplier'; value8 = 0.01; % if you want to increase resolution. CAREFUL!!!!

   
   CS_s = struct(field1,value1,field2,value2,field3,value3,field4,value4,field5,value5,field6,value6,field7,value7,field8,value8);
     



    %% =================== SCAN CURRENT:     ===================
    
   % analysis is the structure which contains the result of the measurement
   % for each loop.
    
   
    [analysis,flag_end_scan] = ...
        scan_current_fromfile_werrors_newDetector_Irene_v3('filein',config_struct.datafile,'cleanfactor',[3,2],'plotmode',0,'loop_index',loop_now,'gui_mode',1,'fig_s',figure_struct,'cs_mode',0, 'CS_struct',CS_s,...
        'Interpol',CS_s.InterpolationMultiplier,'config_s',config_struct); 

    
     %% ===================  Update the meas structure for each loop:  ===================
    
    analysis.loop(1).CS_s = CS_s; % compressed sensing structure
     
    meas.loop(loop_now).analysis      = analysis.loop(1); % raw data
    meas.loop(loop_now).ibase      = ibase_reconstruct(list_meas==1); % ibase
    meas.mean                  = analysis.mean; % averaged data
    meas.loop(loop_now).CS_s          = analysis.loop(1).CS_s; % compressed sensing structure
    meas.loop(loop_now).CS_s.nlist = nlist;
 
    %% =================== COMPRESSED SENSING:     ===================
    % output_struct is the structure which contains the results of the
    % compressed sensing over the whole loop.
   
    index_list_new = CS_s.nlist;
    
    for kk= 1:loop_now 
       % compares the index list of the current sampled with the list for previous loops: 
       index_new = setdiff(meas.loop(kk).CS_s.nlist,index_list_new);       
       index_list_new = sort([index_list_new;index_new]); 
        meas.loop(kk).ibase      = ibase_reconstruct(meas.loop(kk).CS_s.nlist); % raw data
        meas.loop(kk).Preal      = meas.loop(kk).Preal (meas.loop(kk).CS_s.nlist);
        meas.loop(kk).Pimag      = meas.loop(kk).Pimag(meas.loop(kk).CS_s.nlist);
        meas.loop(kk).Pmag      = meas.loop(kk).Pmag(meas.loop(kk).CS_s.nlist);
        meas.loop(kk).deltaPhase = meas.loop(kk).deltaPhase(meas.loop(kk).CS_s.nlist);
    end
    
     %List of the measured currents 

    
     [meas_conc_cs,list_new] = concat_avg_data(meas,loop_now);
     
     CS_s.list = list_new;
     
  %[output_struct] = Pol2ISF(analysis,analysis.mean.Preal(list_new==1),analysis.mean.Pimag(list_new==1),CS_s,list_new);
  [output_struct] = Pol2ISF(meas,meas_conc_cs.Preal,meas_conc_cs.Pimag,CS_s,kk);
  ListHistogram(list_new,config_struct.binvect,CS_s.ibasetoreconstruct,2);

   %% ===================  Update the meas structure for each loop:  ===================
    
   
   analysis.loop(loop_now).CS_out = output_struct; % compressed sensing results
   analysis.list_new = list_new;

    
        
    
    %% =================== Plot results: ===================
    
    % Display histogram of currents
%     bar(handles.current_hist,binranges,bincounts,'histc')
%     axis(handles.current_hist,[1,res_pol,0,1]);
%     title(handles.current_hist,'Sampling Histogram');
%     xlabel(handles.current_hist,'Frequencies');
%     ylabel(handles.current_hist,'Sampling Percentage');


   if min(analysis.ibase)<0                    
        %plot the mean results (temp):
       
        set(figure_struct.plothandle1_1_4,'XData',abs(analysis.ibase(analysis.ibase<0)),'YData',analysis.mean.Preal(analysis.ibase<0));
        set(figure_struct.plothandle1_2_4,'XData',abs(analysis.ibase(analysis.ibase<0)),'YData',analysis.mean.Pimag(analysis.ibase<0));

        set(figure_struct.plothandle1_1_3,'XData',analysis.ibase(analysis.ibase>=0),'YData',analysis.mean.Preal(analysis.ibase>=0));
        set(figure_struct.plothandle1_2_3,'XData',analysis.ibase(analysis.ibase>=0),'YData',analysis.mean.Pimag(analysis.ibase>=0));
        
        
        %set(figure_struct.plothandle4_1_1,'XData',output_struct.time_ps(output_struct.time_ps>0)','YData',real(output_struct.IKt(output_struct.time_ps>0)));
        %set(figure_struct.plothandle4_1_2,'XData',output_struct.time_ps(output_struct.time_ps>0)','YData',imag(output_struct.IKt(output_struct.time_ps>0)));

   else
        set(figure_struct.plothandle1_1_3,'XData',analysis.ibase,'YData',analysis.mean.Preal);
        set(figure_struct.plothandle1_2_3,'XData',analysis.ibase,'YData',analysis.mean.Pimag);
   end


   
   
    %set(figure_struct.plothandle3_1,'XData',output_struct.Energ_meV','YData',real(output_struct.SKw)'/conv_struct.C_J2meV);
    %set(figure_struct.plothandle3_2,'XData',output_struct.Energ_meV','YData',imag(output_struct.SKw)'/conv_struct.C_J2meV);
    
     set(figure_struct.plothandle3_1_2,'XData',NaN,'YData',NaN);
    
    
   set(figure_struct.plothandle3_1,'XData',output_struct.Energ_meV','YData',real(output_struct.SKw)'/output_struct.conv.C_J2meV);
   set(figure_struct.plothandle3_2,'XData',output_struct.Energ_meV','YData',imag(output_struct.SKw)'/output_struct.conv.C_J2meV);
    
    % comparison between the polarization and the calculated IKt
   xdata=analysis.setime;
   ydata_real=analysis.mean.Preal;
   ydata_imag=analysis.mean.Pimag;
   
   XDataPre=output_struct.time_ps(output_struct.time_ps>0)';
   Xlength=length(XDataPre);
   YDataPre=output_struct.IKt(output_struct.time_ps>0);
   
   if CS_s.InterpolationMultiplier>1
        XDataNew=XDataPre(1:floor(Xlength/CS_s.InterpolationMultiplier));
        YDataNew=YDataPre(1:floor(Xlength/CS_s.InterpolationMultiplier));
   else 
        XDataNew=XDataPre;
        YDataNew=YDataPre;
   end
   
   
   set(figure_struct.plothandle4_1_1,'XData',XDataNew*2,'YData',real(YDataNew));
   set(figure_struct.plothandle4_1_2,'XData',XDataNew*2,'YData',imag(YDataNew));
   set(figure_struct.plothandle4_2_1,'XData',xdata(xdata>=0)','YData',ydata_real(xdata>=0));
   set(figure_struct.plothandle4_2_2,'XData',xdata(xdata>=0)','YData',ydata_imag(xdata>=0));
%   
   
%% =================== Grafic part of the GUI: ===================
%  Construct the components for the Preal/Pimag - Pmag/Pphase popup
    %  menu
    hpopup1 = uicontrol('Parent', figure_struct.fh,'Style','popupmenu',...
          'String',{'P_x/P_y','P_{mag}/P_{phase}'},'Position', [330 650 120 20],'Callback',{@popup1_menu_Callback});
      
    %  Construct the components for the log/linear popup menu for
    %  polarization
    hpopup2 = uicontrol('Parent', figure_struct.fh,'Style','popupmenu',...
          'String',{'Linear/Linear','Log/Linear'},'Position', [70 350 120 20],'Callback',{@popup2_menu_Callback});  
      
       %  Construct the components for the log/linear popup menu for
    % IKT
    hpopup3 = uicontrol('Parent', figure_struct.fh,'Style','popupmenu',...
          'String',{'SKw','Linear/Linear','Log/Linear'},'Position', [650 30 120 20],'Callback',{@popup3_menu_Callback}); 
      
    %  Construct the components for the log/linear popup menu for
    % SKw
    hpopup4 = uicontrol('Parent', figure_struct.fh,'Style','popupmenu',...
          'String',{'IKt','Linear/Linear','Log/Linear'},'Position', [70 30 120 20],'Callback',{@popup4_menu_Callback}); 
      
     %  Select to plot the wavelength distribution or the SKw
    hpopup5 = uicontrol('Parent', figure_struct.fh,'Style','popupmenu',...
          'String',{'SKw','\rho(\lambda_F)'},'Position', [670 30 120 20],'Callback',{@popup5_menu_Callback});  

      
    % show figure toolbar:
    set(figure_struct.fh,'toolbar','figure')


   %%  Callbacks for simple_gui. These callbacks automatically
   %  have access to component handles and initialized data 
   %  because they are nested at a lower level.
 
   %  Pop-up menu callback. Read the pop-up menu Value property
   %  to determine which item is currently displayed and make it
   %  the current data.
      function popup1_menu_Callback(source,eventdata) 
         
         % Determine the selected data set.
         str = get(source, 'String');
         val = get(source,'Value');
         
         % Set current data to the selected data set.
         switch str{val};
             case 'P_x/P_y' % User select P_real/P_imag                 
                if min(analysis.ibase)<0                    
                    %plot the mean results (temp):
                   
                    set(figure_struct.plothandle1_1_4,'XData',abs(analysis.ibase(analysis.ibase<0)),'YData',analysis.mean.Preal(analysis.ibase<0));
                    set(figure_struct.plothandle1_2_4,'XData',abs(analysis.ibase(analysis.ibase<0)),'YData',analysis.mean.Pimag(analysis.ibase<0));

                    set(figure_struct.plothandle1_1_3,'XData',analysis.ibase(analysis.ibase>=0),'YData',analysis.mean.Preal(analysis.ibase>=0));
                    set(figure_struct.plothandle1_2_3,'XData',analysis.ibase(analysis.ibase>=0),'YData',analysis.mean.Pimag(analysis.ibase>=0));

                else
                    set(figure_struct.plothandle1_1_3,'XData',analysis.ibase,'YData',analysis.mean.Preal);
                    set(figure_struct.plothandle1_2_3,'XData',analysis.ibase,'YData',analysis.mean.Pimag);
                end

                 xlabel(figure_struct.axhandle1,'Tilted current / A');
                 ylabel(figure_struct.axhandle1,'Real and Imag Polarisation');
                
             case 'P_{mag}/P_{phase}' % User selects P_magnitude/P_phase
                if min(analysis.ibase)<0                    
                    %plot thfigure(10);e mean results (temp):

                    set(figure_struct.plothandle1_1_4,'XData',abs(analysis.ibase(analysis.ibase<0)),'YData',analysis.mean.Pmag(analysis.ibase<0));
                    set(figure_struct.plothandle1_2_4,'XData',abs(analysis.ibase(analysis.ibase<0)),'YData',analysis.mean.deltaPhase(analysis.ibase<0));

                    set(figure_struct.plothandle1_1_3,'XData',analysis.ibase(analysis.ibase>=0),'YData',analysis.mean.Pmag(analysis.ibase>=0));
                    set(figure_struct.plothandle1_2_3,'XData',analysis.ibase(analysis.ibase>=0),'YData',analysis.mean.deltaPhase(analysis.ibase>=0));

                else
                    set(figure_struct.plothandle1_1_3,'XData',analysis.ibase,'YData',analysis.mean.Pmag);
                    set(figure_struct.plothandle1_2_3,'XData',analysis.ibase,'YData',analysis.mean.deltaPhase);
                end
                 
                 xlabel(figure_struct.axhandle1,'Tilted current / A');
                 ylabel(figure_struct.axhandle1,'Modulus and phase of Polarisation');
                
         end
         
      end
  
   function popup2_menu_Callback(source,eventdata) 
         
         % Determine the selected data set.
         str = get(source, 'String');
         val = get(source,'Value');
         
         % Set current data to the selected data set.
         switch str{val};
             case 'Linear/Linear' % User selected linear plot             
                  set(figure_struct.axhandle1,'XScale','Linear')
                 
             case 'Log/Linear' % User selects linear/log plot
                 set(figure_struct.axhandle1,'XScale','Log')
         end
         
   end
  
    function popup3_menu_Callback(source,eventdata) 

             % Determine the selected data set.
             str = get(source, 'String');
             val = get(source,'Value');

             % Set current data to the selected data set.
             switch str{val};
                 case 'Linear/Linear' % User selected to plot the SKW in linear plot             
                      set(figure_struct.axhandle3,'YScale','Linear')

                 case 'Log/Linear' % User selects to plot the SKW linear/log plot
                     set(figure_struct.axhandle3,'YScale','Log')
             end

    end
      
     function popup4_menu_Callback(source,eventdata) 

                 % Determine the selected data set.
                 str = get(source, 'String');
                 val = get(source,'Value');

                 % Set current data to the selected data set.
                 switch str{val};
                     case 'Linear/Linear' % User selected to plot the IKT in linear plot             
                          set(figure_struct.axhandle4,'XScale','Linear')

                     case 'Log/Linear' % User selects to plot the IKT in log/linear plot
                         set(figure_struct.axhandle4,'XScale','Log')
                 end

     end
   
 
 
 
    function popup5_menu_Callback(source,eventdata) 

                 % Determine the selected data set.
                 str = get(source, 'String');
                 val = get(source,'Value');

                 % Set current data to the selected data set.
                 switch str{val};
                     case 'SKw' % User selected to plot the IKT in linear plot             
                          set(figure_struct.plothandle3_1,'XData',output_struct.Energ_meV','YData',real(output_struct.SKw)'/output_struct.conv.C_J2meV);
                          set(figure_struct.plothandle3_2,'XData',output_struct.Energ_meV','YData',imag(output_struct.SKw)'/output_struct.conv.C_J2meV);
                          
                          set(figure_struct.plothandle3_1_2,'XData',NaN,'YData',NaN); % remove interpolated wavelength distribution
                          set(figure_struct.plothandle3_2_2,'XData',NaN,'YData',NaN);
                          
                          xlabel(figure_struct.axhandle3,'\Delta E [meV]');
                          ylabel(figure_struct.axhandle3,'S(\Delta E)_\theta in [1/meV*rad]');
                          title(figure_struct.axhandle3,'Dynamic scattering function at a constant \theta');
                          xlim(figure_struct.axhandle3,[min(output_struct.Energ_meV), 10]); %max(output_struct.Energ_meV)])
                          ylim(figure_struct.axhandle3,[min(real(output_struct.SKw)'/output_struct.conv.C_J2meV), max(real(output_struct.SKw)'/output_struct.conv.C_J2meV)])

                     case '\rho(\lambda_F)' % User selects to plot the IKT in log/linear plot
                        set(figure_struct.plothandle3_1,'XData',output_struct.lambda_m'*1e10,'YData',real(output_struct.Plambda)'*1e-10);
                        set(figure_struct.plothandle3_2,'XData',output_struct.lambda_m'*1e10,'YData',imag(output_struct.Plambda)'*1e-10);
                        set(figure_struct.plothandle3_1_2,'XData',output_struct.lambdaF_m_interp'*1e10,'YData',real(output_struct.Plambda_interp)'*1e-10);
                        set(figure_struct.plothandle3_2_2,'XData',output_struct.lambdaF_m_interp'*1e10,'YData',imag(output_struct.Plambda_interp)'*1e-10);
                        
                        xlabel(figure_struct.axhandle3,'\lambda_F [Angstroms]');
                        ylabel(figure_struct.axhandle3,'\rho(\lambda_F)_\theta in [1/Angstroms*rad]');
                        title(figure_struct.axhandle3,'Final wavelength distribution function at a constant \theta');
                        xlim(figure_struct.axhandle3,[min(output_struct.lambda_m'*1e10) max(output_struct.lambda_m'*1e10)]);
                        ylim(figure_struct.axhandle3,[min(real(output_struct.Plambda)'*1e-10) max(real(output_struct.Plambda)'*1e-10)]);
                        
                 end

    end
   

 
end 