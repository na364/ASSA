classdef Vibrations_cot
    
      %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties        
    end
    
    methods(Static)
        
        function dE = S1(dK)
            x = [0.0000    0.0277    0.0468    0.0702    0.0872    0.1128    0.1340    0.1532    0.1851    0.2021    0.2298    0.2723 ...
                0.3000    0.3319    0.3702    0.4149    0.4575    0.4830    0.5128    0.5490    0.6128    0.6554    0.7129    0.7810 ...
                0.8534    0.8896    0.9407    1.0152    1.0812    1.1174    1.1621    1.2090    1.2665    1.3240    1.3794    1.4071];
            y = [0.0000    0.5078    0.7813    1.0156    1.3672    1.6406    1.9922    2.3828    2.7344    3.0078    3.4375    3.9453 ...
                4.4141    4.8828    5.2344    5.9375    6.2500    6.6406    6.9531    7.2656    8.0078    8.5156    9.1016    9.7266 ...
                10.5078   10.8984   11.2891   11.8359   12.2266   12.6563   12.8125   13.0469   13.1641   13.3984   13.3984   13.3984];
            
            [x,ind] = sort(x); y=y(ind);
            x = [x 2.845-flip(x(1:end-1))];
            y = [y flip(y(1:end-1))];
            dE=interp1(x,y,abs(dK)); 
        end
        
        function dE = LR(dK)
            x = [0.0000    0.0149    0.0340    0.0617    0.0957    0.1340    0.1979    0.2298    0.2851    0.3191 ...
                0.3638    0.4191    0.4553    0.4894    0.5255    0.5702    0.6128    0.6532    0.6915    0.7298 ...
                0.7915    0.8426    0.8894    0.9468    0.9851    1.0362    1.0809    1.1766    1.2660    1.3362    1.4000];
            y = [0.0778    0.6615    1.3230    2.1012    3.1518    4.0078    5.8366    6.7315    8.2490    9.2607 ...
                10.4669   11.7510   12.6070   13.1907   13.8521   14.5914   15.2918   15.9144   16.2646   16.6926 ...
                17.1984   17.7432   18.0545   18.3268   18.4825   18.5992   18.6770   18.7938   18.6381   18.6770 18.6381];
            
            [x,ind] = sort(x); y=y(ind);
            x = [x 2.845-flip(x(1:end-1))];
            y = [y flip(y(1:end-1))];
            dE=interp1(x,y,abs(dK)); 
        end
        
        function v0 = get_v0(Ts,m)
            load_chess_parameters;            
            v0 = sqrt(2*SE_kB*Ts/(m*SE_amu));%m/s
        end
                
        function [fncName,fitFuncStr,fitParams] = get_vibrations(dataSet,E0, Ts, dK, gammaTrue, theta, x, y, refRes)
            
            load_chess_parameters;
            
            if exist('refRes','var') && ~isempty(refRes)
                exp_dcy=refRes.vBest(2)*1000;
                a2_frm_isf = exp_dcy * 1e9 * SE_h / SE_e / pi * 1000;
            end
            
            dE_scancurve = -E0:0.01:20;
            if exist('E0','var') && ~isempty(E0)
                [dk_scancurve]=de2dk(E0,dE_scancurve,gammaTrue,theta);
            end
        
            fitParams = []; %[lower1 upper1 startPoint1; lower2 upper2 startPoint2; ...]
            fncName={};
            fitFuncStr='';
            
            % case 1 - 300K, I/2.55, files dy015085:dy015102, 400K, I/2.55, files dy015103:dy015125
            % case 2 - 300K, I/2, <11(-2)> : files dy0 [14895:14914, 14916:14924, 14979,14981,14983]
            %                     <110> : files dy0 [14941 14943 14945:14956 14958:14964 14966 14967 14969:14973 14975 14976 14978]
            %          200K, I/2, <11(-2)> : files dy0 [15062:15084]
            %          Arrhenius, I/2, <11(-2)> : files dy0 [14925:14939 14985:15008]
            %          Arrhenius, I/2.55, <11(-2)> : files dy0 [15127:15216]
            % case 3 - Simulation
            
            switch dataSet
                case 1
                    
                     % === Jump - Lorenzian
                    fncName(1) = {'cl_Jump'};
                    fitFuncStr = [fitFuncStr '+' '(a1/pi)*((a2/2)/(x^2+(a2/2)^2))'];
                    if dK < 0,  a1=[0 1 0.1]; a2=[0 1 0]; end
                    if exist('a2_frm_isf','var'), a1=[0 1 0.1]; a2=[a2_frm_isf*0.95 a2_frm_isf*1.05 a2_frm_isf]; end
                    fitParams = [fitParams; a1; a2];

%                     fncName(2) = {'cl_Jump2'};
%                     fitFuncStr = [fitFuncStr '+' '(a3/pi)*((a4/2)/(x^2+(a4/2)^2))'];                
%                     if dK < 0,  a3=[0 3 0.7]; a4=[0 1 0]; end
%                     if dK < 0.3,  a3=[0 3 0.7]; a4=[0 1 0]; end
%                     fitParams = [fitParams; a3; a4];


                    % === Balistic - Gaussian ===                    
                    fncName(3) = {'cg_Balistic'};
                    fitFuncStr = [fitFuncStr '+' 'b1*exp(-1*x^2/b2^2)'];
%                     b2_mid = 2*sqrt(log(2))*SE_hbar*abs(dK*1e10)*Vibrations.get_v0(Ts,104.15)/SE_e*1000;
%                     b1 = [0 0.5 0.01]; b2 = [b2_mid*0.95 b2_mid*1.05 b2_mid];
                   if dK<0, b1 = [0 1 0.1]; b2 = [0 10 0.01]; end
                   if dK<-0.225, b1 = [0.1 1 0.1]; b2 = [0.07 10 0.07]; end
                   fitParams = [fitParams; b1; b2];
                   
                   % === Raighly mode - loss and gain Gaussian(s) ===
                   if 0
                       fncName(4)={'ncg_S1'};
                        S1 = @(x) 15.74*sin(0.9272*x)-1.915*sin(0.9272*x).^3;

                        % positive
                        fitFuncStr = [fitFuncStr '+c1*exp(-1*(x-c2)^2/c3^2)'];
                        dE = S1(dk_scancurve);
                        [min_dE, ind]=min(abs(dE-(dE_scancurve)));
                        c2=[dE(ind)-0.7 dE(ind)+0.7 dE(ind)];
                        c1=[0 0 0];
                        if min_dE<1e-1, if -4 < dK && dK < -0.65, c1=[0 0.1 0.05]; end; end
                        if dK<4,  c3 = [0.01 1 0.15]; end
                        fitParams = [fitParams; c1;c2;c3];

                        % negative
                        fitFuncStr = [fitFuncStr '+c4*exp(-1*(x-c5)^2/c6^2)'];
                        dE = -S1(dk_scancurve);
                        [min_dE, ind]=min(abs(dE-(dE_scancurve)));
                        c5=[dE(ind)-0.7 dE(ind)+0.7 dE(ind)];
                        c4=[0 0 0];
                        if min_dE<1e-1, if -4 < dK && dK < -0.65, c4=[0 0.1 0.05]; end; end
                        if dK<4,  c6 = [0.01 1 0.15]; end
                        fitParams = [fitParams; c4;c5;c6];
                   end

                   % === Longitudinal mode - loss and gain Gaussian(s) ===
                   if 0
                        fncName(5)={'ncg_LR'};
                        LR= @(x) 24.12*sin(1.288*x)-5.381*sin(1.288*x).^3;

                        % positive
                        fitFuncStr = [fitFuncStr '+d1*exp(-1*(x-d2)^2/d3^2)'];
                        dE = LR(dk_scancurve);
                        [min_dE, ind]=min(abs(dE-(dE_scancurve)));
                        d2=[dE(ind)-4 dE(ind)+4 dE(ind)];
                        d1=[0 0 0];
                        if min_dE<1e-1, if -4 < dK && dK < -0.65, d1=[0 1 0.5]; end; end
                        if dK<4,  d3 = [0.1 3 0.2]; end
                        fitParams = [fitParams; d1;d2;d3];

                        % negative
                        fitFuncStr = [fitFuncStr '+d4*exp(-1*(x-d5)^2/d6^2)'];
                        dE = -LR(dk_scancurve);
                        [min_dE, ind]=min(abs(dE-(dE_scancurve)));
                        d5=[dE(ind)-4 dE(ind)+4 dE(ind)];
                        d4=[0 0 0];
                        if min_dE<1e-1, if -4 < dK && dK < -0.65, d4=[0 1 0.5]; end; end                        
                        if dK<4,  d6 = [0.1 3 0.2]; end
                        fitParams = [fitParams; d4;d5;d6];
                    end

                    % === isolated mode - loss and gain Gaussian(s) ===
                    if 0
                        fncName(6)={'ncg_isolate1'};                            
                        isolate1 = @(x) repmat(0.65,size(x));

                        % positive
                        fitFuncStr = [fitFuncStr '+e1*exp(-1*(x-e2)^2/e3^2)'];
                        dE = isolate1(dk_scancurve);
                        [min_dE, ind]=min(abs(dE-(dE_scancurve)));
                        e2=[dE(ind)-0.35 dE(ind)+0.35 dE(ind)-0.3];
                        e1=[0 0 0];
                        if min_dE<1e-1, if 0.1<abs(dK) && abs(dK)<0.7, e1=[0 1 0.5]; end; end
                        if dK<4,  e3 = [0.05 0.15 0.05]; end
                        fitParams = [fitParams; e1;e2;e3];

                        % negative
                        fitFuncStr = [fitFuncStr '+e4*exp(-1*(x-e5)^2/e6^2)'];
                        dE = -isolate1(dk_scancurve);
                        [min_dE, ind]=min(abs(dE-(dE_scancurve)));
                        e5=[dE(ind)-0.35 dE(ind)+0.35 dE(ind)+0.3];
                        e4=[0 0 0];
                        if min_dE<1e-1, if 0.1<abs(dK) && abs(dK)<0.7, e4=[0 1 0.5]; end; end
                        if dK<4,  e6 = [0.05 0.15 0.05]; end
                        fitParams = [fitParams; e4;e5;e6];
                    end
                    
                    % === Background - Gaussian(s) ===
                    if 0
                        fncName(7)={'ncg_Background'};                                                
                        fitFuncStr = [fitFuncStr '+f1*exp(-1*(x-f2)^2/f3^2)'];
                        f1 = [0 0 0];
                        if -4 < dK && dK < -3, f1 = [0 0.025 0.01]; end                    
                        f2 = [-0.5 0.5 0]; f3 = [4 20 6];
                        fitParams = [fitParams; f1;f2;f3];
                    end
                    
                    % === Manual interactive - Gaussian(s) ===
                    if 0
                        dE_man=[]; amp_man=[];
                        fh = figure; plot(x,y);
                        [dE_man, amp_man] = ginput(3);                    
                        
                        if length(dE_man)>0
                            fncName(8)={'ncg_manual1'};                                                
                            fitFuncStr = [fitFuncStr '+g1*exp(-1*(x-g2)^2/g3^2)'];
                            g1 = [0 2 0.05];
                            g2 = [dE_man(1)-0.3 dE_man(1)+0.3 dE_man(1)]; g3 = [0 2 0.01];
                            fitParams = [fitParams; g1;g2;g3];
                        end

                        if length(dE_man)>1
                            fncName(9)={'ncg_manual2'};                                                
                            fitFuncStr = [fitFuncStr '+h1*exp(-1*(x-h2)^2/h3^2)'];
                            h1 = [0 4 0.05];                     
                            h2 = [dE_man(2)-2 dE_man(2)+2 dE_man(2)]; h3 = [0 10 1];
                            fitParams = [fitParams; h1;h2;h3];
                        end
                        
                        if length(dE_man)>2
                            fncName(10)={'ncg_manual3'};                                                
                            fitFuncStr = [fitFuncStr '+i1*exp(-1*(x-i2)^2/i3^2)'];
                            i1 = [0 10 0.05];
                            i2 = [0.5 20 dE_man(2)]; i3 = [0 20 1];
                            fitParams = [fitParams; i1;i2;i3];
                        end
                        
                        close(fh)
                    end
                    
                    % === Manual fixed - Gaussian(s) ===
                    if 1      
                            fncName(11)={'ncg_manual4'};
                            fitFuncStr = [fitFuncStr '+j1*exp(-1*(x-j2)^2/j3^2)'];
                            j1 = [0 0 0]; j2 = [-5 5 0.3]; j3 = [0 10 1];
                            if dK<     0  , j1 = [0 0.0 0.0]; j2 = [-5 5 0.3]; j3 = [0 10 1]; end
                            if dK< -0.31  , j1 = [0 0.05 0.0]; j2 = [-5 5 0.3]; j3 = [0 10 1]; end
                            if dK< -0.45  , j1 = [0 0.5 0.01]; j2 = [0.05 5 0.3]; j3 = [0.05 1 0.05]; end
                            if dK< -0.6  , j1 = [0 0 0]; j2 = [0.05 5 0.3]; j3 = [0.05 1 0.05]; end
                            fitParams = [fitParams; j1;j2;j3];
                        
                            fncName(12)={'ncg_manual5'};
                            fitFuncStr = [fitFuncStr '+k1*exp(-1*(x-k2)^2/k3^2)'];    
                            k1 = [0 0 0]; k2 = [-5 5 0.3]; k3 = [0 10 1];
                            if dK< 0      , k1 = [0 0.5 0.03]; k2 = [-5 5 0.3]; k3 = [0 10 1]; end
                            if dK< -0.31  , k1 = [0 0.5 0.03]; k2 = [-5 5 0.3]; k3 = [0 10 1]; end
                            if dK< -0.45  , k1 = [0 0.5 0.03]; k2 = [-5 -0.05 -0.3]; k3 = [0.1 1 0.15]; end
                            fitParams = [fitParams; k1;k2;k3];
                        
                            fncName(13)={'ncg_manual6'};                                                
                            fitFuncStr = [fitFuncStr '+l1*exp(-1*(x-l2)^2/l3^2)'];
                            l1 = [0 0 0]; l2 = [0.1 10 0]; l3 = [0 20 1]; 
                            if dK< -0.45, l1 = [0 0.5 0]; l2 = [-5 -0.1 -1]; l3 = [0 0.05 0.01]; end
                            if dK< -0.61, l1 = [0 0.5 0]; l2 = [4 10 6]; l3 = [0 2 1]; end
                            fitParams = [fitParams; l1;l2;l3];
                            
                            fncName(14)={'ncg_manual7'};                                                
                            fitFuncStr = [fitFuncStr '+m1*exp(-1*(x-m2)^2/m3^2)'];
                            m1 = [0 0 0]; m2 = [0.5 20 4]; m3 = [0 20 1];
                            if dK< -0.61, m1 = [0 0.5 0]; m2 = [4 10 6]; m3 = [0 2 1]; end
                            fitParams = [fitParams; m1;m2;m3];                        
                        
                    end
                    
                    % ===== Background =====
                    fncName(15)={'lin_BG'};
                    fitFuncStr = [fitFuncStr '+y0'];
                    fitParams = [fitParams; 0 max(min(y),0.01) 0];
                
                case 2
                    
                     % === Jump - Lorenzian
                    fncName(1) = {'cl_Jump'};
                    fitFuncStr = [fitFuncStr '+' '(a1/pi)*((a2/2)/(x^2+(a2/2)^2))'];
%                     [exp1] = BProb_baravis(2.55,4,dK,[],112); fwhm_lorentz = (exp1*1e9)*2*SE_hbar*1000/SE_e;
%                     a1=[0 3 0.1]; a2=[fwhm_lorentz*0.9 fwhm_lorentz*1.1 fwhm_lorentz];
                    if dK < 0,  a1=[0 1 0.1]; a2=[0 10 0]; end
                    if exist('a2_frm_isf','var'), a1=[0 1 0.1]; a2=[a2_frm_isf*0.95 a2_frm_isf*1.05 a2_frm_isf]; end
                    fitParams = [fitParams; a1; a2];

%                     fncName(2) = {'cl_Jump2'};
%                     fitFuncStr = [fitFuncStr '+' '(a3/pi)*((a4/2)/(x^2+(a4/2)^2))'];                
%                     if dK < 0,  a3=[0 3 0.7]; a4=[0 1 0]; end
%                     if dK < 0.3,  a3=[0 3 0.7]; a4=[0 1 0]; end
%                     fitParams = [fitParams; a3; a4];


                    % === Balistic - Gaussian ===                    
                    fncName(3) = {'cg_Balistic'};
                    fitFuncStr = [fitFuncStr '+' 'b1*exp(-1*x^2/b2^2)'];
                    fwhm_ballistic = 2*sqrt(log(2))*SE_hbar*abs(dK*1e10)*Vibrations.get_v0(Ts,104.15)/SE_e*1000;
                     b2_mid = fwhm_ballistic/(2*sqrt(2*log(2)));
                    b1 = [0 1 0.01]; b2 = [b2_mid*0 b2_mid*10 b2_mid];
%                    if dK<0, b1 = [0 1 0.1]; b2 = [0 10 0.01]; end
%                    if dK<-0.225, b1 = [0.1 1 0.1]; b2 = [0.07 10 0.07]; end
                    fitParams = [fitParams; b1; b2];

                    % === Raighly mode - loss and gain Gaussian(s) ===
                    if 0
                        fncName(4)={'ncg_S1'};
                        S1 = @Vibrations_cot.S1;

                        % positive
                        fitFuncStr = [fitFuncStr '+c1*exp(-1*(x-c2)^2/c3^2)'];
                        dE = S1(dk_scancurve);
                        [min_dE, ind]=min(abs(dE-(dE_scancurve)));
                        c2=[dE(ind)-0.4 dE(ind)+0.4 dE(ind)];
                        c1=[0 0 0];
                        if min_dE<1e-1, if dK < -0.3, c1=[0 0.5 0.05]; end; end
                        if dK<4,  c3 = [0 2 0.15]; end
                        fitParams = [fitParams; c1;c2;c3];

                        % negative
                        fitFuncStr = [fitFuncStr '+c4*exp(-1*(x-c5)^2/c6^2)'];
                        dE = -S1(dk_scancurve);
                        [min_dE, ind]=min(abs(dE-(dE_scancurve)));
                        c5=[dE(ind)-0.4 dE(ind)+0.4 dE(ind)];
                        c4=[0 0 0];
                        if min_dE<1e-1, if -4 < dK && dK < -0.65, c4=[0 0.5 0.05]; end; end
                        if dK<4,  c6 = [0.01 2 0.15]; end
                        fitParams = [fitParams; c4;c5;c6];
                    end

                    % === Longitudinal mode - loss and gain Gaussian(s) ===
                    if 0
                        fncName(5)={'ncg_LR'};
                        LR= @(x) 24.12*sin(1.288*x)-5.381*sin(1.288*x).^3;

                        % positive
                        fitFuncStr = [fitFuncStr '+d1*exp(-1*(x-d2)^2/d3^2)'];
                        dE = LR(dk_scancurve);
                        [min_dE, ind]=min(abs(dE-(dE_scancurve)));
                        d2=[dE(ind)-4 dE(ind)+4 dE(ind)];
                        d1=[0 0 0];
                        if min_dE<1e-1, if -4 < dK && dK < -0.65, d1=[0 1 0.5]; end; end
                        if dK<4,  d3 = [0.1 3 0.2]; end
                        fitParams = [fitParams; d1;d2;d3];

                        % negative
                        fitFuncStr = [fitFuncStr '+d4*exp(-1*(x-d5)^2/d6^2)'];
                        dE = -LR(dk_scancurve);
                        [min_dE, ind]=min(abs(dE-(dE_scancurve)));
                        d5=[dE(ind)-4 dE(ind)+4 dE(ind)];
                        d4=[0 0 0];
                        if min_dE<1e-1, if -4 < dK && dK < -0.65, d4=[0 1 0.5]; end; end                        
                        if dK<4,  d6 = [0.1 3 0.2]; end
                        fitParams = [fitParams; d4;d5;d6];
                    end

                    % === isolated mode - loss and gain Gaussian(s) ===
                    if 0
                        fncName(6)={'ncg_isolate1'};                            
                        isolate1 = @(x) repmat(0.65,size(x));

                        % positive
                        fitFuncStr = [fitFuncStr '+e1*exp(-1*(x-e2)^2/e3^2)'];
                        dE = isolate1(dk_scancurve);
                        [min_dE, ind]=min(abs(dE-(dE_scancurve)));
                        e2=[dE(ind)-0.35 dE(ind)+0.35 dE(ind)-0.3];
                        e1=[0 0 0];
                        if min_dE<1e-1, if 0.1<abs(dK) && abs(dK)<0.7, e1=[0 1 0.5]; end; end
                        if dK<4,  e3 = [0.05 0.15 0.05]; end
                        fitParams = [fitParams; e1;e2;e3];

                        % negative
                        fitFuncStr = [fitFuncStr '+e4*exp(-1*(x-e5)^2/e6^2)'];
                        dE = -isolate1(dk_scancurve);
                        [min_dE, ind]=min(abs(dE-(dE_scancurve)));
                        e5=[dE(ind)-0.35 dE(ind)+0.35 dE(ind)+0.3];
                        e4=[0 0 0];
                        if min_dE<1e-1, if 0.1<abs(dK) && abs(dK)<0.7, e4=[0 1 0.5]; end; end
                        if dK<4,  e6 = [0.05 0.15 0.05]; end
                        fitParams = [fitParams; e4;e5;e6];
                    end
                    
                    % === Background - Gaussian(s) ===
                    if 0
                        fncName(7)={'ncg_Background'};                                                
                        fitFuncStr = [fitFuncStr '+f1*exp(-1*(x-f2)^2/f3^2)'];
                        f1 = [0 0 0];
                        if -4 < dK && dK < -3, f1 = [0 0.025 0.01]; end                    
                        f2 = [-0.5 0.5 0]; f3 = [4 20 6];
                        fitParams = [fitParams; f1;f2;f3];
                    end
                    
                    % === Manual interactive - Gaussian(s) ===
                    if 0
                        dE_man=[]; amp_man=[];
                        fh = figure; plot(x,y);
                        [dE_man, amp_man] = ginput(3);                    
                        
                        if length(dE_man)>0
                            fncName(8)={'ncg_manual1'};                                                
                            fitFuncStr = [fitFuncStr '+g1*exp(-1*(x-g2)^2/g3^2)'];
                            g1 = [0 2 0.05];
                            g2 = [dE_man(1)-0.3 dE_man(1)+0.3 dE_man(1)]; g3 = [0 2 0.01];
                            fitParams = [fitParams; g1;g2;g3];
                        end

                        if length(dE_man)>1
                            fncName(9)={'ncg_manual2'};                                                
                            fitFuncStr = [fitFuncStr '+h1*exp(-1*(x-h2)^2/h3^2)'];
                            h1 = [0 4 0.05];                     
                            h2 = [dE_man(2)-2 dE_man(2)+2 dE_man(2)]; h3 = [0 10 1];
                            fitParams = [fitParams; h1;h2;h3];
                        end
                        
                        if length(dE_man)>2
                            fncName(10)={'ncg_manual3'};                                                
                            fitFuncStr = [fitFuncStr '+i1*exp(-1*(x-i2)^2/i3^2)'];
                            i1 = [0 10 0.05];
                            i2 = [0.5 20 dE_man(2)]; i3 = [0 20 1];
                            fitParams = [fitParams; i1;i2;i3];
                        end
                        
                        close(fh)
                    end
                    
                    % === Manual fixed - Gaussian(s) ===
                    if 1      
                            fncName_tmp ={'ncg_manual4'};
                            fitFuncStr_tmp =  '+j1*exp(-1*(x-j2)^2/j3^2)';
                            j1 = []; j2 = []; j3 = [];
                            if dK< -0.15 && -0.25 <= dK  , j1 = [0 0.03 0.01]; j2 = [-0.3 -0.2 -0.25]; j3 = [0.05 1 0.05]; fncName(11) = fncName_tmp; fitFuncStr = [fitFuncStr fitFuncStr_tmp]; end
                            if dK< -0.45 && -0.65 <= dK  , j1 = [0 0.5 0.01]; j2 = [0.2 1 0.9]; j3 = [0.05 1 0.05]; fncName(11) = fncName_tmp; fitFuncStr = [fitFuncStr fitFuncStr_tmp]; end
                            if dK< -0.85 && -1.45 <= dK  , j1 = [0 0.5 0.01]; j2 = [2 4.5 4.2]; j3 = [0.05 2 0.05]; fncName(11) = fncName_tmp; fitFuncStr = [fitFuncStr fitFuncStr_tmp]; end
                            if dK< -2.15 && -2.35 <= dK  , j1 = [0 0.5 0.01]; j2 = [0.2 1 0.9]; j3 = [0.05 1 0.05]; fncName(11) = fncName_tmp; fitFuncStr = [fitFuncStr fitFuncStr_tmp]; end
                            if dK< -2.45 && -2.85 <= dK  , j1 = [0 0.5 0.01]; j2 = [0.2 1 0.9]; j3 = [0.05 1 0.05]; fncName(11) = fncName_tmp; fitFuncStr = [fitFuncStr fitFuncStr_tmp]; end
                            fitParams = [fitParams; j1;j2;j3];
                        
                            fncName_tmp={'ncg_manual5'};
                            fitFuncStr_tmp = '+k1*exp(-1*(x-k2)^2/k3^2)';    
                            k1 = []; k2 = []; k3 = [];
                            if dK< -0.25 && -0.45 <= dK , k1 = [0.01 0.5 0.0]; k2 = [-0.9 -0.5 -0.6]; k3 = [0.01 0.3 0.1]; fncName(12) = fncName_tmp; fitFuncStr = [fitFuncStr fitFuncStr_tmp]; end
                            if dK< -0.45 && -0.65 <= dK, k1 = [0 0.5 0.03]; k2 = [-5 -0.05 -0.3]; k3 = [0.1 1 0.15]; fncName(12) = fncName_tmp; fitFuncStr = [fitFuncStr fitFuncStr_tmp]; end
                            fitParams = [fitParams; k1;k2;k3];
                        
                            fncName_tmp={'ncg_manual6'};                                                
                            fitFuncStr_tmp = '+l1*exp(-1*(x-l2)^2/l3^2)';
                            l1 = []; l2 = []; l3 = [];
                            if dK< -0.45 && -0.61 <= dK, l1 = [0 0.5 0]; l2 = [-10 -1 -2]; l3 = [0 5 1]; fncName(13) = fncName_tmp; fitFuncStr = [fitFuncStr fitFuncStr_tmp]; end
                            if dK< -0.61 && -0.71 <= dK, l1 = [0 0.5 0]; l2 = [4 15 4]; l3 = [0 5 1]; fncName(13) = fncName_tmp; fitFuncStr = [fitFuncStr fitFuncStr_tmp]; end
                            if dK< -0.71 && -0.95 <= dK, l1 = [0 0.5 0]; l2 = [5 8 6.4]; l3 = [0 5 1]; fncName(13) = fncName_tmp; fitFuncStr = [fitFuncStr fitFuncStr_tmp]; end
                            if dK< -0.95 && -1.05 <= dK, l1 = [0 0.5 0]; l2 = [7 10 8.8]; l3 = [0 5 1]; fncName(13) = fncName_tmp; fitFuncStr = [fitFuncStr fitFuncStr_tmp]; end
                            if dK< -1.05 && -1.15 <= dK, l1 = [0 0.5 0]; l2 = [9 9.5 9.2]; l3 = [0 5 1]; fncName(13) = fncName_tmp; fitFuncStr = [fitFuncStr fitFuncStr_tmp]; end
                            if dK< -1.15 && -1.25 <= dK, l1 = [0 0.5 0]; l2 = [9.7 10.7 10]; l3 = [0 5 1]; fncName(13) = fncName_tmp; fitFuncStr = [fitFuncStr fitFuncStr_tmp]; end
                            if dK< -1.25 && -1.35 <= dK, l1 = [0.03 0.5 0.1]; l2 = [10.9 11.5 11]; l3 = [0 3 1]; fncName(13) = fncName_tmp; fitFuncStr = [fitFuncStr fitFuncStr_tmp]; end
                            if dK< -1.35 && -1.45 <= dK, l1 = [0 0.5 0]; l2 = [11.5 12.5 12]; l3 = [0 5 1]; fncName(13) = fncName_tmp; fitFuncStr = [fitFuncStr fitFuncStr_tmp]; end
                            if dK< -1.45 && -1.55 <= dK, l1 = [0 0.5 0]; l2 = [12 13.1 12.5]; l3 = [0 5 1]; fncName(13) = fncName_tmp; fitFuncStr = [fitFuncStr fitFuncStr_tmp]; end
                            if dK< -1.55 && -1.65 <= dK, l1 = [0 0.5 0]; l2 = [12.5 14 12.8]; l3 = [0 5 1]; fncName(13) = fncName_tmp; fitFuncStr = [fitFuncStr fitFuncStr_tmp]; end
                            if dK< -1.65 && -2.85 <= dK  , l1 = [0 0.5 0]; l2 = [7 12 9]; l3 = [0 5 1]; fncName(13) = fncName_tmp; fitFuncStr = [fitFuncStr fitFuncStr_tmp]; end
                            fitParams = [fitParams; l1;l2;l3];
                            
                            fncName_tmp={'ncg_manual7'};                                                
                            fitFuncStr_tmp = '+m1*exp(-1*(x-m2)^2/m3^2)';
                            m1 = []; m2 = []; m3 = [];
                            if dK< -0.61 && -0.71 <= dK  , m1 = [0 0.5 0]; m2 = [7 15 8]; m3 = [0 5 1]; fncName(14) = fncName_tmp; fitFuncStr = [fitFuncStr fitFuncStr_tmp]; end
                            if dK< -0.71 && -1.25 <= dK  , m1 = [0 0.5 0]; m2 = [7 15 9.5]; m3 = [0 5 1]; fncName(14) = fncName_tmp; fitFuncStr = [fitFuncStr fitFuncStr_tmp]; end
                            if dK< -1.25 && -1.35 <= dK  , m1 = [0 0.5 0]; m2 = [12.5 18 13]; m3 = [0 10 1]; fncName(14) = fncName_tmp; fitFuncStr = [fitFuncStr fitFuncStr_tmp]; end
                            if dK< -1.35 && -2.85 <= dK  , m1 = [0 0.5 0]; m2 = [7 15 9.5]; m3 = [0 5 1]; fncName(14) = fncName_tmp; fitFuncStr = [fitFuncStr fitFuncStr_tmp]; end
                            if dK< -2.85 && -3.05 <= dK  , m1 = [0 0.5 0]; m2 = [3 15 4]; m3 = [0 5 1]; fncName(14) = fncName_tmp; fitFuncStr = [fitFuncStr fitFuncStr_tmp]; end
                            fitParams = [fitParams; m1;m2;m3];                        
                        
                    end
                    
                    % ===== Background =====
                    fncName(15)={'lin_BG'};
                    fitFuncStr = [fitFuncStr '+y0'];
                    fitParams = [fitParams; 0 max(min(y),0.01) 0];
                
                case 3 % for simulation
                    
                     % === Jump - Lorenzian
                    fncName(1) = {'cl_Jump'};
                    fitFuncStr = [fitFuncStr '+' '(a1/pi)*((a2/2)/(x^2+(a2/2)^2))'];
%                     [ISF1, exp1] = BProb_baravis(2.55,4,0.1,0,112); exp1 = exp1*2*SE_hbar;
%                     a1=[0 3 0.1]; a2=[exp1*0.95 exp1*1.05 exp1];
                    if dK >= 0,  a1=[0 2 0.1]; a2=[0 10 0.1]; end
                    if exist('a2_frm_isf','var'), a1=[0 2 0.1]; a2=[a2_frm_isf*0.95 a2_frm_isf*1.05 a2_frm_isf]; end
                    fitParams = [fitParams; a1; a2];

%                     fncName(2) = {'cl_Jump2'};
%                     fitFuncStr = [fitFuncStr '+' '(a3/pi)*((a4/2)/(x^2+(a4/2)^2))'];                
%                     if dK < 0,  a3=[0 3 0.7]; a4=[0 1 0]; end
%                     if dK >= 0,  a3=[0 3 0.7]; a4=[0 1 0]; end
%                     fitParams = [fitParams; a3; a4];


%                     % === Balistic - Gaussian ===                    
                    fncName(3) = {'cg_Balistic'};
                    fitFuncStr = [fitFuncStr '+' 'b1*exp(-1*x^2/b2^2)'];                    
                    fwhm_ballistic = 2*sqrt(log(2))*SE_hbar*abs(dK*1e10)*Vibrations.get_v0(Ts,104.15)/SE_e*1000;
                    b2_mid = fwhm_ballistic/(2*sqrt(log(2)));
                    b1 = [0 1 0.5]; b2 = [b2_mid*0.99 b2_mid*1.01 b2_mid*1];
%                    if dK>=0, b1 = [0 1 0.1]; b2 = [0 10 0.01]; end
                    fitParams = [fitParams; b1; b2];

                    % === Raighly mode - loss and gain Gaussian(s) ===
                    if 0
                        fncName(4)={'ncg_S1'};
                        S1 = @(x) 15.74*sin(1.288*x)-1.915*sin(1.288*x).^3;

                        % positive
                        fitFuncStr = [fitFuncStr '+c1*exp(-1*(x-c2)^2/c3^2)'];
                        dE = S1(dk_scancurve);
                        [min_dE, ind]=min(abs(dE-(dE_scancurve)));
                        c2=[dE(ind)-0.7 dE(ind)+0.7 dE(ind)];
                        c1=[0 0 0];
                        if min_dE<1e-1, if -4 < dK && dK < -0.65, c1=[0 0.1 0.05]; end; end
                        if dK<4,  c3 = [0.01 1 0.15]; end
                        fitParams = [fitParams; c1;c2;c3];

                        % negative
                        fitFuncStr = [fitFuncStr '+c4*exp(-1*(x-c5)^2/c6^2)'];
                        dE = -S1(dk_scancurve);
                        [min_dE, ind]=min(abs(dE-(dE_scancurve)));
                        c5=[dE(ind)-0.7 dE(ind)+0.7 dE(ind)];
                        c4=[0 0 0];
                        if min_dE<1e-1, if -4 < dK && dK < -0.65, c4=[0 0.1 0.05]; end; end
                        if dK<4,  c6 = [0.01 1 0.15]; end
                        fitParams = [fitParams; c4;c5;c6];
                    end

                    % === Longitudinal mode - loss and gain Gaussian(s) ===
                    if 0
                        fncName(5)={'ncg_LR'};
                        LR= @(x) 24.12*sin(1.288*x)-5.381*sin(1.288*x).^3;

                        % positive
                        fitFuncStr = [fitFuncStr '+d1*exp(-1*(x-d2)^2/d3^2)'];
                        dE = LR(dk_scancurve);
                        [min_dE, ind]=min(abs(dE-(dE_scancurve)));
                        d2=[dE(ind)-4 dE(ind)+4 dE(ind)];
                        d1=[0 0 0];
                        if min_dE<1e-1, if -4 < dK && dK < -0.65, d1=[0 1 0.5]; end; end
                        if dK<4,  d3 = [0.1 3 0.2]; end
                        fitParams = [fitParams; d1;d2;d3];

                        % negative
                        fitFuncStr = [fitFuncStr '+d4*exp(-1*(x-d5)^2/d6^2)'];
                        dE = -LR(dk_scancurve);
                        [min_dE, ind]=min(abs(dE-(dE_scancurve)));
                        d5=[dE(ind)-4 dE(ind)+4 dE(ind)];
                        d4=[0 0 0];
                        if min_dE<1e-1, if -4 < dK && dK < -0.65, d4=[0 1 0.5]; end; end                        
                        if dK<4,  d6 = [0.1 3 0.2]; end
                        fitParams = [fitParams; d4;d5;d6];
                    end

                    % === isolated mode - loss and gain Gaussian(s) ===
                    if 0
                        fncName(6)={'ncg_isolate1'};                            
                        isolate1 = @(x) repmat(0.65,size(x));

                        % positive
                        fitFuncStr = [fitFuncStr '+e1*exp(-1*(x-e2)^2/e3^2)'];
                        dE = isolate1(dk_scancurve);
                        [min_dE, ind]=min(abs(dE-(dE_scancurve)));
                        e2=[dE(ind)-0.35 dE(ind)+0.35 dE(ind)-0.3];
                        e1=[0 0 0];
                        if min_dE<1e-1, if 0.1<abs(dK) && abs(dK)<0.7, e1=[0 1 0.5]; end; end
                        if dK<4,  e3 = [0.05 0.15 0.05]; end
                        fitParams = [fitParams; e1;e2;e3];

                        % negative
                        fitFuncStr = [fitFuncStr '+e4*exp(-1*(x-e5)^2/e6^2)'];
                        dE = -isolate1(dk_scancurve);
                        [min_dE, ind]=min(abs(dE-(dE_scancurve)));
                        e5=[dE(ind)-0.35 dE(ind)+0.35 dE(ind)+0.3];
                        e4=[0 0 0];
                        if min_dE<1e-1, if 0.1<abs(dK) && abs(dK)<0.7, e4=[0 1 0.5]; end; end
                        if dK<4,  e6 = [0.05 0.15 0.05]; end
                        fitParams = [fitParams; e4;e5;e6];
                    end
                    
                    % === Background - Gaussian(s) ===
                    if 0
                        fncName(7)={'lin_Background'};                                                
                        fitFuncStr = [fitFuncStr '+f1*exp(-1*(x-f2)^2/f3^2)'];
                        f1 = [0 1 0];
                        if -4 < dK && dK < -3, f1 = [0 0.025 0.01]; end                    
                        f2 = [-0.5 0.5 0]; f3 = [0 20 1];
                        fitParams = [fitParams; f1;f2;f3];
                    end
                    
                    % === Manual interactive - Gaussian(s) ===
                    if 0
                        dE_man=[]; amp_man=[];
                        fh = figure; plot(x,y);
                        [dE_man, amp_man] = ginput(3);                    
                        
                        if length(dE_man)>0
                            fncName(8)={'ncg_manual1'};                                                
                            fitFuncStr = [fitFuncStr '+g1*exp(-1*(x-g2)^2/g3^2)'];
                            g1 = [0 2 0.05];
                            g2 = [dE_man(1)-0.3 dE_man(1)+0.3 dE_man(1)]; g3 = [0 2 0.01];
                            fitParams = [fitParams; g1;g2;g3];
                        end

                        if length(dE_man)>1
                            fncName(9)={'ncg_manual2'};                                                
                            fitFuncStr = [fitFuncStr '+h1*exp(-1*(x-h2)^2/h3^2)'];
                            h1 = [0 4 0.05];                     
                            h2 = [dE_man(2)-2 dE_man(2)+2 dE_man(2)]; h3 = [0 10 1];
                            fitParams = [fitParams; h1;h2;h3];
                        end
                        
                        if length(dE_man)>2
                            fncName(10)={'ncg_manual3'};                                                
                            fitFuncStr = [fitFuncStr '+i1*exp(-1*(x-i2)^2/i3^2)'];
                            i1 = [0 10 0.05];
                            i2 = [0.5 20 dE_man(2)]; i3 = [0 20 1];
                            fitParams = [fitParams; i1;i2;i3];
                        end
                        
                        close(fh)
                    end
                    
                    % === Manual fixed - Gaussian(s) ===
                    if 1      
%                             fncName(11)={'ncg_manual4'};
%                             fitFuncStr = [fitFuncStr '+j1*exp(-1*(x-j2)^2/j3^2)'];
%                             j1 = [0 1 0.01]; j2 = [-5 -2.5 -3]; j3 = [0 1 0.05];
%                             fitParams = [fitParams; j1;j2;j3];
%                         
%                             fncName(12)={'ncg_manual5'};
%                             fitFuncStr = [fitFuncStr '+k1*exp(-1*(x-k2)^2/k3^2)'];    
%                             k1 = [0 1 0.01]; k2 = [2.5 5 3]; k3 = [0 1 0.05];               
%                             fitParams = [fitParams; k1;k2;k3];
%                         
%                             fncName(13)={'ncg_manual6'};                                                
%                             fitFuncStr = [fitFuncStr '+l1*exp(-1*(x-l2)^2/l3^2)'];
%                             l1 = [0 1 0.01]; l2 = [-8 -3 -6]; l3 = [0 1 0.05];
%                             fitParams = [fitParams; l1;l2;l3];
%                             
%                             fncName(14)={'ncg_manual7'};                                                
%                             fitFuncStr = [fitFuncStr '+m1*exp(-1*(x-m2)^2/m3^2)'];
%                             m1 = [0 1 0.01]; m2 = [3 8 6]; m3 = [0 1 0.05];
%                             fitParams = [fitParams; m1;m2;m3];
                            
                              % T-mode with NON uniform width, and some slack in the harmonics.
                              % The higher harmonics have with which
                              % depends on the lower harmonic. Same with
                              % the frequencies
                            fncName(15) = {'sp1_manual8'};
                            fitFuncStr = [fitFuncStr '+' '(n1/pi)*         ((n2/2)/         ((x-1*n3)^2+      (n2/2)^2))'];
                            fitFuncStr = [fitFuncStr '+' '(n1/pi)*         ((n2/2)/         ((x+1*n3)^2+      (n2/2)^2))'];
                            fitFuncStr = [fitFuncStr '+' '(n4*n1/pi)*      ((n5*n2/2)/      ((x-(2*n3+n10))^2+(n5*n2/2)^2))'];
                            fitFuncStr = [fitFuncStr '+' '(n4*n1/pi)*      ((n5*n2/2)/      ((x+(2*n3+n10))^2+(n5*n2/2)^2))'];
                            fitFuncStr = [fitFuncStr '+' '(n6*n4*n1/pi)*   ((n7*n5*n2/2)/   ((x-(3*n3+n11))^2+(n7*n5*n2/2)^2))'];
                            fitFuncStr = [fitFuncStr '+' '(n6*n4*n1/pi)*   ((n7*n5*n2/2)/   ((x+(3*n3+n11))^2+(n7*n5*n2/2)^2))'];
                            fitFuncStr = [fitFuncStr '+' '(n8*n6*n4*n1/pi)*((n9*n7*n5*n2/2)/((x-(4*n3+n12))^2+(n9*n7*n5*n2/2)^2))'];
                            fitFuncStr = [fitFuncStr '+' '(n8*n6*n4*n1/pi)*((n9*n7*n5*n2/2)/((x+(4*n3+n12))^2+(n9*n7*n5*n2/2)^2))'];
                            if dK >= 0
                                n1=[0 3 0]; n4=[0.001 0.5 0.1]; n6=[0.001 0.5 0.1]; n8=[0.001 0.5 0.1];
                                n2=[0 10 1]; n5=[1 5 1];n7=[1 5 1]; n9=[1 5 1];
                                n3=[-2.2 -1.95 -2]; n10=[-2 0 0]; n11=[-2 0 0]; n12=[-2 0 0];
                            end
                            fitParams = [fitParams; n1; n10; n11; n12; n2; n3; n4; n5; n6; n7; n8; n9];
                            
                            % T-mode with fixed width, and some slack in the harmonics.
%                             fncName(16) = {'sp2_manual9'};
%                             fitFuncStr = [fitFuncStr '+' '(n1/pi)*((n2/2)/((x-1*n3)^2+(n2/2)^2))'];
%                             fitFuncStr = [fitFuncStr '+' '(n1/pi)*((n2/2)/((x+1*n3)^2+(n2/2)^2))'];
%                             fitFuncStr = [fitFuncStr '+' '(n4*n1/pi)*((n2/2)/((x-(2*n3+n7))^2+(n2/2)^2))']; % 2nd harmonic (positive)
%                             fitFuncStr = [fitFuncStr '+' '(n4*n1/pi)*((n2/2)/((x+(2*n3+n7))^2+(n2/2)^2))']; % 2nd harmonic (negative)
%                             fitFuncStr = [fitFuncStr '+' '(n5*n4*n1/pi)*((n2/2)/((x-(3*n3+n8))^2+(n2/2)^2))']; % 3rd harmonic (positive)
%                             fitFuncStr = [fitFuncStr '+' '(n5*n4*n1/pi)*((n2/2)/((x+(3*n3+n8))^2+(n2/2)^2))']; % 3rd harmonic (negative)
%                             fitFuncStr = [fitFuncStr '+' '(n6*n5*n4*n1/pi)*((n2/2)/((x-(4*n3+n9))^2+(n2/2)^2))']; % 4th harmonic (positive)
%                             fitFuncStr = [fitFuncStr '+' '(n6*n5*n4*n1/pi)*((n2/2)/((x+(4*n3+n9))^2+(n2/2)^2))']; % 4th harmonic (negative)
%                             if dK >= 0
%                                 n1=[0.001 3 0.1]; n2=[0 10 1]; n3=[-2.2 -2.05 -2.1];
%                                 n4=[0.04 0.5 0.1]; n5=[0.04 0.5 0.1];
%                                 n6=[0.04 0.5 0.1];
%                                 n7=[-0.5 0 0]; % allow minor shift (since present in the simulation)
%                                 n8=[-0.5 0 0]; % allow minor shift (since present in the simulation)
%                                 n9=[-0.5 0 0]; % allow minor shift (since present in the simulation)                            
%                             end
%                             fitParams = [fitParams; n1; n2; n3; n4; n5; n6; n7; n8; n9];

                    end
                    
                    % ===== Background =====
                    fncName(17)={'lin_BG'};
                    fitFuncStr = [fitFuncStr '+y0'];
                    fitParams = [fitParams; 0 max(min(y),0.01) 0];
                otherwise
            end

        end
        
    end
end