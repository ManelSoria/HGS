%***********************************************************************************************************
% *HGS 2.0
% *By Caleb Fuster, Manel Soria and Arnau Miró
% *ESEIAAT UPC               
%***********************************************************************************************************
%
% Adiabatic H2O2 decomposition
%
% Inlet: H2O2
% Outlet: H2O + (1/2)O2 at Tp

clear; 
 
Tr=300;         % K
p=1;            % not relevant
 
% Enthalpy of each species 
hH2O2=@(T) HGSsingle('H2O2','h',T,p);
hO2=@(T)  HGSsingle('O2','h',T,p);
hH2O=@(T) HGSsingle('H2O','h',T,p);

% Equation to be solved
eq=@(Tp) hH2O2(Tr)-...
         (hH2O(Tp)+0.5*hO2(Tp));
 
options=optimset(...
        'Display','none',...
        'MaxIter',1000,...
        'TolFun', 1.0e-4);
 
 
[Tp,fval,exitflag]=fzero(eq,3000,options);

fprintf('Temperature of the products: %.2f K \n',Tp);

