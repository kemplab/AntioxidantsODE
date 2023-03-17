function time_series = blap_ODE_system(endtime, step, celldens, protein_abund, blap_conc, parameter_pctchange, species_pctchange)
    % model for compartmentalized hydrogen peroxide consumption
    
    
    % Establish genes of interest and index
    genes = ["GPX1", "GPX2", "GPX3", "GPX4", "GPX5", "GPX6", "GPX7", "GPX8", "PRDX1", "PRDX2", ...
        "PRDX3", "PRDX6", "CAT", "TXN", "TXN2", "TXNRD1", "TXNRD2", "TXNRD3", "GLRX", "GLRX2", ...
        "G6PD", "GLUD1", "GSR", "GSTP1", "POR", "NQO1", "SOD1", "SOD2", "SOD3", "AQP3", ...
        "AQP8", "AQP9", "GCLC"];

    % Solver Parameters
    ti = 0; % start time
    tf = endtime; % stop time  
    tspan = [ti tf];
    options=odeset('AbsTol',1E-10,'RelTol',1E-4,'maxstep',step); %,'Jacobian',@J
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Defining forward reaction rates ()
        % Defining the array, k, that will house reaction rates
        k = zeros(36,1);
 
 %      Rxns (Cytosolic)
    
        % H2O2 Transport
        % Permeability constant (plasma membrane)
        % H2O2cyt -> H2O2out 
        k(1) = 2e-4*protein_abund(30)/0.0302; % cm/s
        
        % Mitochondrial production of H2O2in
        k(2) = 4;
        
        % GPx1_red reacting with H2O2in
        % H2O2cyt + GPX1r -> GPX1o ; EC 1.11.1.9
        % Updated 9/22 with Kassi Stein PLoS Comp Bio paper k2
        k(3) = 60; % uM-1*s-1 (60 uM-1*s-1 Kassi paper)
        
        % GPx1_ox reacting with GSH
        % GPX1o + GSH -> GPX1-SG
        % Updated 9/22 with KS PLoS Comp Bio paper k3
        k(4) = 0.04; % uM-1*s-1 
        
        % GPx1-SG reacting with GSH
        % GPX1-SG + GSH -> GPX1r + GSSG
        % Updated 9/22 with KS PLoS Comp Bio Paper k4
        k(5) = 10; % uM-1*s-1 
        
        % Catalase reacts with H2O2 inside peroxisome
        % H2O2p -> 
        % Updated 10/6 
        k(6) = 34; % uM-1*s-1    
        
        % Km of NADP+ 
        % Updated 9/22 with KS PLoS Comp Bio Paper k4
        k(7) = 57; % uM 
        
        % Peroxiredoxin is oxidized by H2O2 
        % UPdated 10/6
        k(8) = 40; % uM-1*s-1
        
        % Oxidized peroxiredoxin is over-oxidized by H2O2 
        % Updated 10/6
        k(9) = 0.072; % uM-1*s-1  
        
        % Reduction of overoxidized Prx by Srx enzyme
        k(10) = 3e-3; % s-1
        
        % Self-catalyzed disulfide formation of Prx-SS from Prx-SOH
        k(11) = 15; % s-1  
        
        % Peroxiredoxin is reduced by Thioredoxin 
        % Updated 10/6
        k(12) = 2.1; % uM-1*s-1 
        
        % Auto-oxidation of GSH 
        k(13) = 7.4e-05; % s-1 
        
        % Protein monothiol oxidized by H2O2
        % Updated 10/6
        k(14) = 0.01; % uM-1*s-1 
        
        % Oxidized Protein monothiol glutathionylated by GSH
        % Updated 10/6
        k(15) = 0.12; % uM-1*s-1
        
        % Grx-SH de-glutathionylates Protein-SSG
        % Updated 10/6
        k(16) = 0.091; % uM-1*s-1 
        
        % GSH de-glutathionylates Grx-SSG
        % Updated 10/6        
        k(17) = 0.037; % uM-1*s-1 
       
        % Protein dithiol oxidized by H2O2
        % Updated 10/6
        k(18) = 5e-1; % uM-1*s-1 
           
        % Protein disulfide reduced by Trx1 
        % Updated 10/6
        k(19) = 1e-1; % uM-1*s-1 
        
        % GSSG reduced by GR
        k(20) = 3.2*protein_abund(23)/0.1419; % uM-1*s-1 
        
        % Oxidized Thioredoxin reduced by TrxR  
        k(21) = 2e1*protein_abund(16)/0.1816; % uM-1*s-1 
        
        % Production of NADPH by G6P-DH
        k(22) = 3.75e2*(protein_abund(21)+protein_abund(22))/(0.2581+0.5159); % uM/s     
        
        % Permeability constant (peroxisomal membrane)
        k(23) = 3e-3; % cm/s 
        
        % GSH synthesis
        k(24) = 4.1e-1*protein_abund(33)/.0075; % uM/s 4.1e-7 
        
        % GSSG efflux 
        k(25) = 1.2e-2; % uM/s 1.2e-8
        
        % GSH + GSSG efflux 
        k(26) = 1.2e-1; % uM/s 1.2e-7
        
        % Trx_SH efflux 
        %k(27) = 7.45e-4; % uM/s 7.45e-10
        
        % Trx_SH synthesis
        %k(28) = 6.97e-4; % uM/s 6.97e-10
        
        % NQO1 Reaction; b-lap Q->HQ
        k(29) = 161.500000*protein_abund(26)/.2905; % 1/uM-s
        
        % b-lap HQ->SQ
        k(30) = 3e2; % 1/uM-s
        
        % b-lap SQ->Q
        k(31) = 3e2; % 1/uM-s
        
        % SOD
        k(32) = 6400.000000*protein_abund(27)/1.6641; %1/uM-s
        
        % CPR
        k(33) = 1900*protein_abund(25)/0.2231; % 1/uM-s   
        
        % blap perm
        k(34) = 1e-6; % cm/s
        
        % blapHQ + GSH -> blapHQ-SG
        k(35) = 0.05*protein_abund(24)/2.6052; % 1/uM-s
        
        % blapHQ-SG -> blapHQ-SG extracell
        k(36) = 1e-6; % cm/s
        
        k= k.*parameter_pctchange;
            
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
        % Defining initial Conditions (Concentrations in uM = micromoles/Liter)
        % Defining the array, x0, that will house initial concentrations
        x0 = zeros(32,1);

        %%% Cytoplasm
        
        % H2O2media
        x0(1) = 0;
        
        % H2O2cytosol
        x0(2) = 1.0e-3; 
       
        % GPXr 
        x0(3) = protein_abund(1);
        
        % GPXo
        x0(4) = x0(3)/5e9;
        
        % GPX-SG
        x0(5) = x0(3)/5e9;      
        
        % GSH
        x0(6) = 368; 
        
        % GSSG
        x0(7) = x0(6)/206;   
        
        % Catalase
        x0(8) = protein_abund(13);  
        
        % H2O2peroxisome
        x0(9) = 1.0e-10;   
        
        % Prx1/2-SH
        x0(10) = protein_abund(9) + protein_abund(10);
        
        % Prx1/2-SOH
        x0(11) = x0(10)/1.92e9;   
        
        % Prx1/2-SOOH
        x0(12) = x0(10)/1.92e9; 
        
        % Prx1/2-SS 
        x0(13) = x0(10)*(.5/100);     
        
        % Trx1-SH
        x0(14) = protein_abund(14);
        
        % Trx1-SS
        x0(15) = x0(14)/5.7;    
        
        % Pr-SH
        x0(16) = 50;%1.22e2; 
        
        % Pr-SOH
        x0(17) = x0(16)*(.5/100);
        
        % Pr-SSG
        x0(18) = x0(16)*(.5/100);  
      
        % Grx-SH 
        x0(19) = protein_abund(19); 
        
        % Grx-SSG
        x0(20) = x0(19)*(.5/100);   
        
        % Pr-(SH)2 
        x0(21) = 450;  
        
        % Pr-SS
        x0(22) = x0(21)*(.5/100);    
        
        % NADPH 
        x0(23) = 3.0e1; 
        
        % NADP+
        x0(24) = x0(23)/100;   
        
        % blap Q
        x0(25) = 0;
        
        % blap HQ
        x0(26) = 0;
        
        % blap SQ
        x0(27) = 0;
        
        % superO2
        x0(28) = 0;
        
        % oxygen
        x0(29) = 2.6e2; 
        
        % blap extracell
        x0(30) = blap_conc;
        
        % blap-GSH
        x0(31) = 0;
        
        % blap-GSH extracell
        x0(32) = 0;
        
        x0 = x0.*species_pctchange;
        
        x0(4) = x0(3)/5e9;     % GPXo 
        x0(5) = x0(3)/5e9;     % GPX-SG 
        x0(7) = x0(6)/206;  % GSSG 
        x0(11) = x0(10)/1.92e9;               % Prx-SOH 
        x0(12) = x0(10)/1.92e9;               % Prx-SOOH 
        x0(13) = x0(10)*(.5/100);     % Prx-SS 
        x0(15) = x0(14)/5.7;   % Trx1-SS 
        x0(17) = x0(16)*(.5/100);  % Pr-SOH
        x0(18) = x0(16)*(.5/100);  % Pr-SSG
        x0(20) = x0(19)*(.5/100);  % Grx-SSG 
        x0(22) = x0(21)*(.5/100);    % Pr-SS
        x0(24) = x0(23)/100;  % NADP+
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % integration
        %tic
        [t, x]=ode15s(@ode_sys,tspan,x0,options,k, celldens);
        %toc
        time_series = [t, x];
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % description of derivatives

        function dxdt=ode_sys(~, x, k, celldens)
                           
             dxdt= x;  % set up vector containing derivatives
                 
             x_area = celldens*1*1000*4*pi()*.00001^2; % in cm^2
             c_area = celldens*1*4*pi()*.001^2 - x_area; % in cm^2
             x_vol = celldens*1*1000*4/3*pi()*.00001^3; % in cm^3
             c_vol = celldens*1*4/3*pi()*.001^3 - x_vol; % in cm^3
             m_vol = 1000 - x_vol - c_vol; % in cm^3

             %Cytosolic Dynamics
                    
             dxdt(1) = k(1)*c_area*(x(2)- x(1))/m_vol; % H2O2 extracellular

             dxdt(2) = k(1)*c_area*(x(1)-x(2))/c_vol - k(23)*x_area*(x(2)-x(9))/c_vol + k(2) - k(3)*x(3)*(x(2)- x0(2))  - k(8)*x(10)*(x(2)- x0(2)) - k(9)*x(11)*(x(2)- x0(2)) - k(14)*x(16)*(x(2)- x0(2)) - k(18)*x(21)*(x(2)- x0(2)) + k(32)*x(28)^2; % H2O2 cytosol    

             dxdt(3) = - k(3)*x(3)*(x(2)- x0(2)) + k(5)*(x(5)-x0(5))*x(6); % GPXr 

             dxdt(4) = k(3)*x(3)*(x(2)- x0(2)) - k(4)*(x(4)-x0(4))*x(6); % GPXo 

             dxdt(5) = k(4)*(x(4)-x0(4))*x(6) - k(5)*(x(5)-x0(5))*x(6); % GPX-SG     

             dxdt(6) = -k(4)*(x(4)-x0(4))*x(6) - k(5)*(x(5)-x0(5))*x(6) - 2*k(13)*x(6) - k(15)*(x(17)-x0(17))*x(6) - k(35)*x(26)*x(6) - k(17)*(x(20)-x0(20))*x(6) + 2*k(20)*(x(7)-x0(7))*x(23) + k(24) - k(26)*x(6)*x(7); % GSH

             dxdt(7) = k(5)*(x(5)-x0(5))*x(6) + k(13)*x(6) + k(17)*(x(20)-x0(20))*x(6) - k(20)*(x(7)-x0(7))*x(23) - k(25) - k(26)*x(6)*x(7);  % GSSG 

             dxdt(8) = 0; % FeCat  

             dxdt(9) = k(23)*x_area*(x(2)-x(9))/x_vol - k(6)*x(8)*(x(9)-x0(9)); % H2O2 peroxisome  

             dxdt(10) = - k(8)*x(10)*(x(2)- x0(2)) + k(12)*(x(13)-x0(13))*x(14); % Prx-SH    

             dxdt(11) = k(8)*x(10)*(x(2)- x0(2)) - k(9)*x(11)*(x(2)- x0(2)) + k(10)*(x(12)-x0(12)) - k(11)*x(11); % Prx-SOH    

             dxdt(12) = k(9)*x(11)*(x(2)- x0(2)) - k(10)*(x(12)-x0(12)); % Prx-SOOH  

             dxdt(13) = k(11)*x(11) - k(12)*(x(13)-x0(13))*x(14); % Prx-SS  

             dxdt(14) = - k(12)*(x(13)-x0(13))*x(14) - k(19)*(x(22)-x0(22))*x(14) + k(21)*(x(15)-x0(15))*x(23); % Trx-SH   

             dxdt(15) = k(12)*(x(13)-x0(13))*x(14) + k(19)*(x(22)-x0(22))*x(14) - k(21)*(x(15)-x0(15))*x(23); % Trx-SS   

             dxdt(16) = - k(14)*x(16)*(x(2)- x0(2)) + k(16)*x(19)*(x(18)-x0(18)); % Pr-SH  

             dxdt(17) = k(14)*x(16)*(x(2)- x0(2)) - k(15)*(x(17)-x0(17))*x(6); % Pr-SOH  

             dxdt(18) = k(15)*(x(17)-x0(17))*x(6) - k(16)*x(19)*(x(18)-x0(18)); % Pr-SSG  

             dxdt(19) = k(17)*(x(20)-x0(20))*x(6) - k(16)*x(19)*(x(18)-x0(18)); % Grx-SH   

             dxdt(20) = k(16)*x(19)*(x(18)-x0(18)) - k(17)*(x(20)-x0(20))*x(6); % Grx-SSG   

             dxdt(21) = - k(18)*x(21)*(x(2)- x0(2)) + k(19)*(x(22)-x0(22))*x(14); % Pr-(SH)2   

             dxdt(22) = k(18)*x(21)*(x(2)- x0(2)) - k(19)*(x(22)-x0(22))*x(14); % Pr-SS   

             dxdt(23) = - k(20)*(x(7)-x0(7))*x(23) - k(21)*(x(15)-x0(15))*x(23) + k(22)*(x(24)-x0(24))/(k(7) + x(24)) - k(29)*x(25)*x(23); % NADPH    

             dxdt(24) = k(20)*(x(7)-x0(7))*x(23) + k(21)*(x(15)-x0(15))*x(23) - k(22)*(x(24)-x0(24))/(k(7) + x(24)) + k(29)*x(25)*x(23); % NADP+

             dxdt(25) = k(34)*c_area*(x(30)- x(25))/c_vol + k(31)*x(27)*x(29) - k(29)*x(25)*x(23) - k(33)*x(25)*x(23);% b-lap Q cyto

             dxdt(26) = k(29)*x(25)*x(23) - k(30)*x(26)*x(29) - k(35)*x(26)*x(6) + k(33)*x(27)*x(23); % b-lap HQ

             dxdt(27) = k(30)*x(26)*x(29) - k(31)*x(27)*x(29) + k(33)*x(25)*x(23) - k(33)*x(27)*x(23);% b-lap SQ

             dxdt(28) = k(30)*x(26)*x(29) + k(31)*x(27)*x(29) - k(32)*x(28)^2;% superO2

             dxdt(29) = - k(30)*x(26)*x(29) - k(31)*x(27)*x(29) + k(32)*x(28)^2;%O2

             dxdt(30) = -k(34)*(x(30)- x(25))/m_vol; % b-lap Q extracellular

             dxdt(31) = k(35)*x(26)*x(6) - k(36)*c_area*(x(31)-x(32))/c_vol; % b-lapHQ-GSH
             
             dxdt(32) = k(36)*c_area*(x(31)-x(32))/m_vol; % b-lapHQ-GSH extracellular
                
        end

end