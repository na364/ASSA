classdef Vibrations_GRu
    
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods(Static)
        
        function dE = RW(dK)
            x = [0, 0.034, 0.079, 0.121, 0.173, 0.223, 0.276, 0.346, 0.414, 0.493, 0.582, 1.0, 1.5, 2];
            y = [0, 0.872, 1.738, 2.732, 3.670, 4.886, 6.139, 7.291, 8.573, 9.849, 11.04, 17, 19, 19];
            
            dE=interp1(x,y,abs(dK));
        end
        
        function dE = LR(dK)
            x = [0;0.0883561643835614;0.248630136986301;...
                0.400684931506849;0.517808219178082;0.573287671232876];
            y = [0;5.26561924001357;11.7124688499934;...
                15.4362844124630;18.1550201276966;19.3961027470988];
            
            dE=interp1(x,y,abs(dK));
        end
        
        function dE = S7(dK)
            x = [-0.00445523869215489;0.0494933667606810;0.112347226574742;...
                0.155029771231589;0.204486576830670;0.267256868275358;...
                0.368202235453881;0.433207980779275;0.509584247362373;...
                0.549958215815314;0.843988300428288];
            y = [0.128886799679650;3.78664821198510;6.99388035795118;...
                9.56031721160207;12.9614540896271;15.1406560116996;...
                19.4369145861625;21.6159249973885;26.1712977471360;...
                27.8383996657265;42.3984905463282];
            
            dE=interp1(x,y,abs(dK));
        end
        
        function dE = S6(dK)
            x = [0.00225634597304921;0.0588843622688811;0.146960200564086;...
                0.205771440509767;0.271289042097566;0.357181656742923;...
                0.402303353180821];
            y = [0.192564156133562;9.31154462202722;20.2911487168773;...
                28.7674187819910;37.2431143145653;48.8654288101953;...
                53.9374978237404];
            
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
            
            % case 1 - 300K, Graphene on Ru(0001)
            
            switch dataSet
                case 1
                    
                    % === Balistic - Gaussian ===
                    if 1
                        fncName(1) = {'cl_Jump'};
                        fitFuncStr = [fitFuncStr '+' '(a1/pi)*((a2/2)/(x^2+(a2/2)^2))'];
%                         b2_mid = 2*sqrt(log(2))*SE_hbar*abs(dK*1e10)*Vibrations.get_v0(Ts,104.15)/SE_e*1000;
%                         b1 = [0 0.5 0.01]; b2 = [b2_mid*0.95 b2_mid*1.05 b2_mid];
                        a1 = [0 1 0.1]; a2 = [0 10 0.01];
                        fitParams = [fitParams; a1; a2];
                    end
                    
                    % === Raighly mode - loss and gain Gaussian(s) ===
                    if 1
                        fncName(2)={'ncg_{RW}'};
                        RW = @Vibrations_GRu.RW;
                        
                        % positive
                        fitFuncStr = [fitFuncStr '+b1*exp(-1*(x-b2)^2/b3^2)'];
                        dE = RW(dk_scancurve);
                        [~, ind]=min(abs(dE-(dE_scancurve)));
                        b2=[dE(ind)-1 dE(ind)+3 dE(ind)*(1+0.16)];
                        b1=[0 0.1 0.005];
                        b3=[0 2 1];
                        if dK>-0.35&&dK<-0.25, b3=[0,1.5 0.5]; end
                        fitParams = [fitParams; b1;b2;b3];
                    end
                    
                    %negative RW mode
                    if 0
                        fncName(3)={'ncg_{nRW}'};
                        RW = @Vibrations_GRu.RW;
                        % negative
                        fitFuncStr = [fitFuncStr '+b4*exp(-1*(x-b5)^2/b6^2)'];
                        dE = -RW(dk_scancurve);
                        [~, ind]=min(abs(dE-(dE_scancurve)));
                        b5=[dE(ind)-0.7 dE(ind)+0.7 dE(ind)];
                        b4=[0 0.1 0.005];
                        b6=[0 2 1];
                        fitParams = [fitParams; b4;b5;b6];
                    end
                    
                    if 0
                        %the peak on the right of RW mode
                        fncName(4)={'ncg_{rofRW}'};
                        RW = @Vibrations_GRu.RW;
                        fitFuncStr = [fitFuncStr '+b7*exp(-1*(x-b8)^2/b9^2)'];
                        dE = RW(dk_scancurve);
                        [~, ind]=min(abs(dE-(dE_scancurve)));
                        b8=[dE(ind) dE(ind)*(1+0.35) dE(ind)*(1+0.2)];
                        b7=[0 0.1 0.005];
                        b9=[0 2 1];
                        fitParams = [fitParams; b7;b8;b9];
                    end
                    
                    if 1
                        fncName(5)={'ncg_{LR}'};
                        LR = @Vibrations_GRu.LR;
                        
                        % positive
                        fitFuncStr = [fitFuncStr '+c1*exp(-1*(x-c2)^2/c3^2)'];
                        dE = LR(dk_scancurve);
                        [~, ind]=min(abs(dE-(dE_scancurve)));
                        c2=[dE(ind)-3 dE(ind)+3 dE(ind)];
                        c1=[0 0.1 0.005];
                        c3=[0 9 1.5];
                        if dK>-0.15, c2=[3 5 4]; end
                        if dK<-0.35&&dK>-0.45, c2=[3 5 4.3]; end
                        fitParams = [fitParams; c1;c2;c3];
                    end
                        
                    if 0
                        fncName(6)={'ncg_{nLR}'};
                        LR = @Vibrations_GRu.LR;
                        % negative
                        fitFuncStr = [fitFuncStr '+c4*exp(-1*(x-c5)^2/c6^2)'];
                        dE = -LR(dk_scancurve);
                        [~,ind]=min(abs(dE-(dE_scancurve)));figure;plot(dk_scancurve,dE);hold on; plot(dk_scancurve,dE_scancurve);
                        c5=[dE(ind)-1.7 dE(ind)+1.7 dE(ind)];disp(dE(ind));
                        c4=[0 0.1 0.001];
                        c6=[0 1 0.15];
                        fitParams = [fitParams; c4;c5;c6];
                    end
                    
                    % === Longitudinal mode - loss and gain Gaussian(s) ===
                    if dK>-0.85&&dK<-0.25
                        fncName(7)={'ncg_{S7}'};
                        S7= @Vibrations_GRu.S7;
                        
                        % positive
                        fitFuncStr = [fitFuncStr '+d1*exp(-1*(x-d2)^2/d3^2)'];
                        dE = S7(dk_scancurve);
                        [~, ind]=min(abs(dE-(dE_scancurve)));
                        d2=[dE(ind)-2 dE(ind)+2 dE(ind)];
                        d1=[0 0.1 0.001];
                        d3 = [0 5 2];
                        fitParams = [fitParams; d1;d2;d3];
                    end
                    if 0
                        % negative S7 mode
                        fncName(8)={'ncg_{nS7}'};
                        S7= @Vibrations_GRu.S7;
                        
                        fitFuncStr = [fitFuncStr '+d4*exp(-1*(x-d5)^2/d6^2)'];
                        dE = -S7(dk_scancurve);
                        [~, ind]=min(abs(dE-(dE_scancurve)));
                        d5=[dE(ind)-2 dE(ind)+2 dE(ind)];
                        d4=[0 0.1 0.001];
                        d6=[0 3 0.2];
                        fitParams = [fitParams; d4;d5;d6];
                    end
                    
                    if dK>-0.45&&dK<-0.25
                        fncName(9)={'ncg_{S6}'};
                        S6= @Vibrations_GRu.S6;
                        fitFuncStr = [fitFuncStr '+e1*exp(-1*(x-e2)^2/e3^2)'];
                        dE = S6(dk_scancurve);
                        [~, ind]=min(abs(dE-(dE_scancurve)));
                        e2=[dE(ind)-2 dE(ind)+2 dE(ind)];
                        if dK<-0.65, e2=[dE(ind)*(1-0.3) dE(ind)*(1+0.3) dE(ind)]; end
                        e1=[0 1 0.001];
                        e3=[0 6 2];
                        fitParams = [fitParams; e1;e2;e3];
                    end
                                        
                    if 0
                        fncName(10)={'ncg_{nS6}'};
                        S6= @Vibrations_GRu.S6;
                        fitFuncStr = [fitFuncStr '+e4*exp(-1*(x-e5)^2/e6^2)'];
                        dE = -S6(dk_scancurve);
                        [~, ind]=min(abs(dE-(dE_scancurve)));
                        e5=[dE(ind)-0.7 dE(ind)+0.7 dE(ind)];
                        e4=[0 0.1 0.001];
                        e6=[0 6 2];
                        fitParams = [fitParams; e4;e5;e6];
                    end
                                        
                    if 1
                        fncName(11)={'ncg_{frtoEla}'};%far right to the elastic
                        fitFuncStr = [fitFuncStr '+f1*exp(-1*(x-f2)^2/f3^2)'];
                        f1 = [0 0.1 0.001];
                        f2 = [0 8 2];
                        f3 = [0 20 1];
                        if dK>-0.15, f2=[6 12 10]; end
                        if dK<-0.65, f2=[3,9,6.5]; end
                        if dK<-0.15&&dK>-0.25, f1=[0 0.001 0.0005]; f2=[2 4 3]; end
                        fitParams = [fitParams; f1;f2;f3];
                    end
                    
                    if dK<-0.4001||dK>-0.3995
                        fncName(12)={'ncg_{nrtoEla}'};%near right to the elastic
                        fitFuncStr = [fitFuncStr '+g1*exp(-1*(x-g2)^2/g3^2)'];
                        g1 = [0 0.1 0.002];
                        g2 = [-1 2 0.5];
                        g3 = [0 2 0.2];
                        if dK>-0.15, g2=[6 10 8.2];g3 = [0 2 1]; end
                        if dK<-0.15&&dK>-0.35, g2=[-1 0 -0.5]; g3=[0 2 0.3]; end
                        if dK<-0.35&&dK>-0.4, g1 = [0.001 0.1 0.003];g2 = [0 2 1.3]; g3 = [0.01 2 0.2]; end
                        fitParams = [fitParams; g1;g2;g3];
                    end
                    
                    % === Background - Gaussian(s) ===
                    if dK<-0.35
                        fncName(13)={'ncg_{bg}'};
                        fitFuncStr = [fitFuncStr '+h1*exp(-1*(x-h2)^2/h3^2)'];
                        h1 = [0 0.1 0.001];
                        h2 = [10 50 30]; 
                        h3 = [0 40 16];
                        fitParams = [fitParams; h1;h2;h3];
                    end
                    
                    % === peaks on dE<0 (negative peaks) - Gaussian(s) ===
                    if 1
                        fncName(14)={'ncg_{np1}'};
                        fitFuncStr = [fitFuncStr '+i1*exp(-1*(x-i2)^2/i3^2)'];
                        i1 = [0 0.1 0.05];
                        i2 = [-5 -0.5 -2.8]; 
                        i3 = [0 10 0.3];
                        if dK<-0.15&&dK>-0.25, i1 = [0 0.1 0.0003]; i2=[-2 -0.2 -1]; i3 = [0 10 0.7];end
                        if dK<-0.25&&dK>-0.45, i2=[-3 -0.2 -1.3]; end
                        if dK<-0.55, i2=[-5 -0.2 -1.4]; i3=[0 10 1]; end
                        
                        fitParams = [fitParams; i1;i2;i3];
                    end
                    
                    % === peaks on dE<0 (negative peaks) - Gaussian(s) ===
                    if dK<-0.75||(dK>-0.65&&dK<-0.4001)||(dK>-0.3995&&dK<-0.35)||dK>-0.15
                        fncName(15)={'ncg_{np2}'};
                        fitFuncStr = [fitFuncStr '+j1*exp(-1*(x-j2)^2/j3^2)'];
                        j1 = [0 0.1 0.05];
                        j2 = [-5 -0.5 -2.5];
                        j3 = [0 10 0.3];
                        if dK>-0.15, j2=[-3 -1 -1.7]; j3 = [0 10 0.9]; end
                        if dK<-0.35&&dK>-0.45, j2=[-3 -1 -2.2]; j3=[0 10 0.4];end
                        if dK<-0.55, j2=[-5 -2 -4.5]; j3=[0 10 1];end
                        fitParams = [fitParams; j1;j2;j3];
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
                    if 0
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
                    fncName(16)={'lin_BG'};
                    fitFuncStr = [fitFuncStr '+y0'];
                    if dK>-0.15
                        fitParams = [fitParams; 0 max(min(y),0.00002) 0];
                    else
                        fitParams = [fitParams; 0 max(min(y),0.002) 0];
                    end
                case 2
                    
                    % === Balistic - Gaussian ===
                    if 1
                        fncName(1) = {'cl_{Jump}'};
                        fitFuncStr = [fitFuncStr '+' '(a1/pi)*((a2/2)/(x^2+(a2/2)^2))'];
%                         b2_mid = 2*sqrt(log(2))*SE_hbar*abs(dK*1e10)*Vibrations.get_v0(Ts,104.15)/SE_e*1000;
%                         b1 = [0 0.5 0.01]; b2 = [b2_mid*0.95 b2_mid*1.05 b2_mid];
                        a1 = [0 1 0.1]; a2 = [0 10 0.01];
                        fitParams = [fitParams; a1; a2];
                    end
                    
                    % === Raighly mode - loss and gain Gaussian(s) ===
                    if 1
                        fncName(2)={'ncg_{RW}'};
                        RW = @Vibrations_GRu.RW;
                        
                        % positive
                        fitFuncStr = [fitFuncStr '+b1*exp(-1*(x-b2)^2/b3^2)'];
                        dE = RW(dk_scancurve);
                        [~, ind]=min(abs(dE-(dE_scancurve)));
                        b2=[dE(ind)*0.95 dE(ind)*1.05 dE(ind)];
                        b1=[0 0.1 0.05];
                        b3=[0 6 1];
                        if dK>-0.15, b2=[0 2 0.9]; end
                        if dK>-0.65, b3=[0 1 0.4]; end
                        fitParams = [fitParams; b1;b2;b3];
                    end
                    
                    if 1
                        fncName(3)={'ncg_{postRW}'};
                        RW = @Vibrations_GRu.RW;
                        
                        % positive
                        fitFuncStr = [fitFuncStr '+e1*exp(-1*(x-e2)^2/e3^2)'];
                        dE = RW(dk_scancurve);
                        [~, ind]=min(abs(dE-(dE_scancurve)));
                        e2=[dE(ind) dE(ind)*1.3 dE(ind)*1.1];
                        e1=[0 0.1 0.0005];
                        e3=[0.05 2 0.5];
                        if dK>-0.15, e2=[0 5 4]; end
                        if dK<-0.75, e1=[0 1 0.005]; end
                        fitParams = [fitParams; e1;e2;e3];
                    end
                    
                    if 1
                        fncName(4)={'ncg_{main}'};
                        RW = @Vibrations_GRu.RW;
                        
                        % positive
                        fitFuncStr = [fitFuncStr '+d1*exp(-1*(x-d2)^2/d3^2)'];
                        dE = RW(dk_scancurve);
                        [~, ind]=min(abs(dE-(dE_scancurve)));
                        d2=[0 10 5];
                        d1=[0.0002 0.1 0.0005];
                        d3=[0.05 15 5];
                        if dK>-0.15, d2=[0 6 4]; end
                        if dK<-0.65, d1=[0 1 0.005]; d3=[0.05 15 5]; end
%                         if dK<-0.85, d3=[0 40 20]; end
                        fitParams = [fitParams; d1;d2;d3];
                    end
                    
                    % === Manual fixed - Gaussian(s) ===
                    if 1
                        fncName(10)={'ncg_{manual1}'};
                        fitFuncStr = [fitFuncStr '+c1*exp(-1*(x-c2)^2/c3^2)'];
                        c1 = [0 1 0.1]; c2 = [0 40 10]; c3 = [0 50 5];
                        if dK<-0.45, c2=[0 40 20]; end
                        if dK<-0.85, c2=[0 50 30]; end
                        fitParams = [fitParams; c1;c2;c3];
                    end
                    
                    if 1
                        fncName(11)={'ncg_{manual2}'};
                        fitFuncStr = [fitFuncStr '+j1*exp(-1*(x-j2)^2/j3^2)'];
                        j1 = [0 1 0.2]; j2 = [-5 0 -1.5]; j3 = [0 50 3];
                        
                        fitParams = [fitParams; j1;j2;j3];
                    end
                    % ===== Background =====
                    fncName(16)={'lin_BG'};
                    fitFuncStr = [fitFuncStr '+y0'];
                    
                    fitParams = [fitParams; 0 max(min(y),0.002) 0];
                case 3 % use one Lorentzian as the RW mode, one Gaussian as the LR mode,
                    % and one Gaussian as the multiphonon background
                    
                    % === Raighly mode - loss and gain Gaussian(s) ===
                    if 1
                        fncName(2)={'ncl_{RW}'};
                        % positive
                        fitFuncStr = [fitFuncStr '+(a1/pi)*((a3/2)/((x-a2)^2+(a3/2)^2))'];
                        a2=[0,20,3.76];
                        a1=[0 0.1 0.05];
                        a3=[0 6 0.15];
                        fitParams = [fitParams; a1;a2;a3];
                        
                    end
                    if 0
                        fncName(2)={'ncg_{RW}'};
                        % positive
                        fitFuncStr = [fitFuncStr '+a1*exp(-1*(x-a2)^2/a3^2)'];
                        b2=[0 10 3.75];
                        b1=[0 0.1 0.01];
                        b3=[0 14 0.2];
                        fitParams = [fitParams; b1;b2;b3];
                    end
                    
                    if 0
                        fncName(3)={'ncl_{postRW}'};
                        % positive
                        fitFuncStr = [fitFuncStr '+(c1/pi)*((c3/2)/((x-c2)^2+(c3/2)^2))'];
                        c2=[0,6,4];%c2=[4.0,4.2,4.0];
                        c1=[0 0.1 0.02];
                        c3=[0 2 0.1];
                        fitParams = [fitParams; c1;c2;c3];
                    end
                    
                    if 1
                        fncName(3)={'ncg_{LR}'};
                        % positive
                        fitFuncStr = [fitFuncStr '+b1*exp(-1*(x-b2)^2/b3^2)'];
                        b2=[4.4 4.6 4.5];
                        b1=[0 0.1 0.01];
                        b3=[0 14 1];
                        fitParams = [fitParams; b1;b2;b3];
                    end
                   
                    
                    if 0
                        fncName(4)={'ncl_{LR}'};
                        % positive
                        fitFuncStr = [fitFuncStr '+(c1/pi)*((c3/2)/((x-c2)^2+(c3/2)^2))'];
                        c2=[4.4 4.6 4.5];
                        c1=[0 0.1 0.01];
                        c3=[0 15 1];
                        fitParams = [fitParams; c1;c2;c3];
                    end
                                        
                    if 1
                        fncName(5)={'ncg_{mph}'};
                        fitFuncStr = [fitFuncStr '+d1*exp(-1*(x-d2)^2/d3^2)'];
                        d2=[0 15 6];
                        d1=[0 0.1 0.001];
                        d3=[0 15 1];
                        fitParams = [fitParams; d1;d2;d3];
                    end
                    if 0
                        fncName(8)={'pl4_{main}'};
                        % a fourth order polynomial
                        fitFuncStr = [fitFuncStr '+n5+n1*x+n2*x.^2+n3*x.^3+n4*x.^4'];
                        n2=[-Inf Inf 0];
                        n1=[-Inf Inf 0];
                        n3=[-Inf Inf 0];
                        n4=[-Inf Inf 0];
                        n5=[-Inf Inf 0];
                        fitParams = [fitParams; n1;n2;n3;n4;n5];
                    end
                    
                    % ===== Background =====
%                     fncName(16)={'lin_BG'};
%                     fitFuncStr = [fitFuncStr '+y0'];
%                     
%                     fitParams = [fitParams; 0 0.1 0];
                case 4 % use one Lorentzian as the RW mode, one Gaussian as the LR mode,
                    % and two Gaussians as the multiphonon background
                    
                    % === Raighly mode - loss and gain Gaussian(s) ===
                    if 1
                        fncName(2)={'ncl_{RW}'};
                        % positive
                        fitFuncStr = [fitFuncStr '+(a1/pi)*((a3/2)/((x-a2)^2+(a3/2)^2))'];
                        a2=[0,20,3.76];
                        a1=[0 0.1 0.05];
                        a3=[0 6 0.15];
                        fitParams = [fitParams; a1;a2;a3];
                        
                    end
                    if 0
                        fncName(2)={'ncg_{RW}'};
                        % positive
                        fitFuncStr = [fitFuncStr '+a1*exp(-1*(x-a2)^2/a3^2)'];
                        b2=[0 10 3.75];
                        b1=[0 0.1 0.01];
                        b3=[0 14 0.2];
                        fitParams = [fitParams; b1;b2;b3];
                    end
                    
                    if 0
                        fncName(3)={'ncl_{postRW}'};
                        % positive
                        fitFuncStr = [fitFuncStr '+(c1/pi)*((c3/2)/((x-c2)^2+(c3/2)^2))'];
                        c2=[0,6,4];%c2=[4.0,4.2,4.0];
                        c1=[0 0.1 0.02];
                        c3=[0 2 0.1];
                        fitParams = [fitParams; c1;c2;c3];
                    end
                    
                    if 1
                        fncName(3)={'ncg_{LR}'};
                        % positive
                        fitFuncStr = [fitFuncStr '+b1*exp(-1*(x-b2)^2/b3^2)'];
                        b2=[4.4 4.6 4.5];
                        b1=[0 0.1 0.01];
                        b3=[0 14 1];
                        fitParams = [fitParams; b1;b2;b3];
                    end
                   
                    
                    if 1
                        fncName(4)={'ncg_{mph1}'};
                        % positive
                        fitFuncStr = [fitFuncStr '+c1*exp(-1*(x-c2)^2/c3^2)'];
                        c2=[3 4.6 4];
                        c1=[0 0.1 0.001];
                        c3=[0 15 1];
                        fitParams = [fitParams; c1;c2;c3];
                    end
                                        
                    if 1
                        fncName(5)={'ncg_{mph2}'};
                        fitFuncStr = [fitFuncStr '+d1*exp(-1*(x-d2)^2/d3^2)'];
                        d2=[0 15 6];
                        d1=[0 0.1 0.001];
                        d3=[0 15 1];
                        fitParams = [fitParams; d1;d2;d3];
                    end
                    
                    % ===== Background =====
%                     fncName(16)={'lin_BG'};
%                     fitFuncStr = [fitFuncStr '+y0'];
%                     
%                     fitParams = [fitParams; 0 0.1 0];

                case 5 % use one Lorentzian as the RW mode, one Gaussian as the LR mode,
                    % and a polynomial as the multiphonon background
                    
                    % === Raighly mode - loss and gain Gaussian(s) ===
                    if 1
                        fncName(2)={'ncl_{RW}'};
                        % positive
                        fitFuncStr = [fitFuncStr '+(a1/pi)*((a3/2)/((x-a2)^2+(a3/2)^2))'];
                        a2=[0,20,3.76];
                        a1=[0 0.1 0.05];
                        a3=[0 6 0.15];
                        fitParams = [fitParams; a1;a2;a3];
                        
                    end
                    
                    if 1
                        fncName(3)={'ncg_{LR}'};
                        % positive
                        fitFuncStr = [fitFuncStr '+b1*exp(-1*(x-b2)^2/b3^2)'];
                        b2=[4.4 4.6 4.5];
                        b1=[0 0.1 0.01];
                        b3=[0 14 1];
                        fitParams = [fitParams; b1;b2;b3];
                    end
                                        
                    if 1
                        fncName(5)={'pl4_{mph}'};
                        fitFuncStr = [fitFuncStr '+d5+d1*x+d2*x.^2+d3*x.^3+d4*x.^4;'];
                        
                        d1=[-15 15 0];
                        d2=[-15 15 0];
                        d3=[-15 15 0];
                        d4=[-15 15 0];
                        d5=[-15 15 0];
                        fitParams = [fitParams; d1;d2;d3;d4;d5];
                    end
            end
        end
    end
end