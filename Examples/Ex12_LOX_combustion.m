%***********************************************************************************************************
% *HGS 2.0
% *By Caleb Fuster, Manel Soria and Arnau Miró
% *ESEIAAT UPC              
%***********************************************************************************************************
%
% LOX - LH2 combustion
%
% Temperature of stoichiometric combustion of H2 O2
% assuming that the reactives enter inlet the chamber as satured 
% liquid at 10 bar
%
% This example requieres INIST


clear; 

species={'H2','O2' , 'H2O','H','O','OH'};
nr=[2;1;0;0;0;0]'; % mol
P=10;

% In INIST we find the enthalpy of O2 at the inlet:
hO2_INIST=INIST('O2','hl_p',10); % kJ/kg
% and convert it to kJ/mol
hO2_INIST=hO2_INIST*HGSsingle('O2','Mm')/1000; % kJ/mol

% However, this enthalpy is not in the same reference as in HGS
% BUT, the difference of enthalpies is roughly the same (the difference is
% due to the fluid model).
% We choose an arbitrary reference state (10 bar, 300K), accesible to both codes
hO2_ref_INIST=INIST('O2','h_pt',10,300) * HGSsingle('O2','Mm')/1000; % kJ/mol

deltaH_O2=hO2_ref_INIST - hO2_INIST; 

hO2_ref_HGS=HGSsingle('O2','h',300,10);

hO2_inlet_HGS = hO2_ref_HGS - deltaH_O2; % kJ/mol

% We repeat the method for H2
hH2_INIST=INIST('H2','hl_p',10)*HGSsingle('H2','Mm')/1000; % kJ/mol
hH2_ref_INIST=INIST('H2','h_pt',10,300) * HGSsingle('H2','Mm')/1000; % kJ/mol
deltaH_H2=hH2_ref_INIST - hH2_INIST; 
hH2_ref_HGS=HGSsingle('H2','h',300,10);
hH2_inlet_HGS = hH2_ref_HGS - deltaH_H2; % kJ/mol

% Finally, we evaluate the TOTAL inlet enthalpy in the HGS reference:
HinLIQ= 2*hH2_inlet_HGS+1*hO2_inlet_HGS;


% HGStp call, with known H at the inlet
[Tpl,~,~,~] = HGStp(species,nr,'H',HinLIQ,P);

% HGStp call, with known T at the inlet (gas)
[Tpg,~,~,~] = HGStp(species,nr,'T',300,P);


% Printing results
fprintf('HGStp Tpl=%f saturated liquid at 10bar inlet K \n',Tpl);
fprintf('HGStp Tpg=%f 300K gas inlet K\n',Tpg);

