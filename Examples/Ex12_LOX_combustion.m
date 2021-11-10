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
% Requires INIST this example
% O2 (NIST) hv(404.36 K)-hl(119.62 K)=14.3753 kJ/mol
% H2 (NIST) hv(413.96K)-hl(31.39K)=10.9495 kJ/mol

function [Tp,np]=Ex12_LOX_combustion

clear; clc;

species={'H2','O2' , 'H2O','H','O','OH'};
nr=[2;1;0;0;0;0]'; % mol
P=10;

% Enthalpy of liquid O2 at Tsat 10 bar (kJ/mol)
h1_O2 = INIST('O2','h_pt',10,119.62);
href_O2 = INIST('O2','h_pt',10,404.36);
hO2= HGSsingle('O2','h',404.36,10) - (href_O2-h1_O2)/(1000*HGSsingle('O2','Mm')); 

% Enthalpy of liquid H2 at Tsat 10 bar (kJ/mol)
h1_H2 = INIST('H2','h_pt',10,31.39);
href_H2 = INIST('H2','h_pt',10,413.96);
hH2= HGSsingle('H2','h',413.96,10) - (href_H2-h1_H2)/(1000*HGSsingle('H2','Mm')); 

% Enthalpy of stoichiometric mixture
HinLIQ=2*hH2+1*hO2;


% HGStp call
[Tp1,n,species,~] = HGStp(species,nr,'H',HinLIQ,P);


% That is what HGStp does
opt = optimoptions(@fsolve,'Display','iter','FunctionTolerance',1e-8,'OptimalityTolerance',1e-8);
Tp=fsolve(@(T)DeltaH(T)-HinLIQ,3000,opt);      % Solver, in our cas we use a HGSsecant
[~,np,~]= HGSeq(species,nr,Tp,P); % After getting Tp, mols equilibrium in a T P

    function DeltaH=DeltaH(Tprod)
        T = Tprod;
        % The result must be in thermodynamic equilbrium at the T suposed
        [~,comp,~] = HGSeq(species,nr,T,P); 
        
        % Enthalpy of the combustion at the T suposed
        Hout = HGSprop(species,comp,T,P,'H');
        
        % Has to be 0 for solving. Pre combustion and Post combusiton
        % enthalpy must be equal
        DeltaH = (Hout);
    end


% Printing results
fprintf('   HGStp <%f>   DeltaH<%f>\n',Tp1,Tp)
for ii=1:length(n)
    fprintf('%s       <%f>             <%f>\n',species{ii},n(ii),np(ii))
end

end