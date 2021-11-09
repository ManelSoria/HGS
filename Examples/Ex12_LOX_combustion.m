%***********************************************************************************************************
% *HGS 2.0
% *By Caleb Fuster, Manel Soria and Arnau Miró
% *ESEIAAT UPC              
%***********************************************************************************************************
%
% LOX - LH2 combustion
%
% Temperature of stoichiometric combustion of H2 O2
% reactives inlet as satured liquit at 10 bar
%
% O2 (NIST) hv(404.36 K)-hl(119.62 K)=14.3753 kJ/mol
% H2 (NIST) hv(413.96K)-hl(31.39K)=10.9495 kJ/mol

function [Tp,np]=Ex12_LOX_combustion

clear;

species={'H2','O2' , 'H2O','H','O','OH'};
nr=[2;1;0;0;0;0]; % mol
P=10;


options = [];

[Tgas,ngas,~,~]= HGStp(species,nr,'T',300,P,options);


% Enthalpy of liquid O2 at Tsat 10 bar (kJ/mol)
hO2= HGSsingle('O2','h',404.36,10) -14.3753; 

% Enthalpy of liquid H2 at Tsat 10 bar (kJ/mol)
hH2= HGSsingle('H2','h',413.96,10) -10.9495; 

% Enthalpy of stoichiometric mixture
HinLIQ=2*hH2+1*hO2;

[Tp1,n,species,flag] = HGStp(species,nr,'H',HinLIQ,P);

 
    function DeltaH=DeltaH(Tprod)
        [~,comp,~] = HGSeq(species,nr,Tprod,P);
        Hout = HGSprop(species,comp,Tprod,P,'H');
        DeltaH = Hout-HinLIQ;
    end

opt = optimset('Display','iter'); 
Tp=fzero(@DeltaH,3000,opt);
[~,np,~]= HGSeq(species,nr,Tp,P);

fprintf('Inlets are gas at 300K: Tp= %f \n',Tgas)
fprintf('Inlets are saturated liquid at 10bar: Tp= %f \n',Tp1)
for ii=1:length(n)
    fprintf('%s \t %f \t %f\n',species{ii},ngas(ii),np(ii))
end

end