%***********************************************************************************************************
%* HGS 2.0 (adapted from the original HGS 1.3) 
%* Original by Arnau Miro, Pau Manent, Manel Soria 
%* Adapted by Caleb Fuster
%
%* ESEIAAT UPC                
%***********************************************************************************************************
%
% Example 02: Adiabatic H2O2 decomposition
%
% Inlet: H2O2
% Outlet: H2O + (1/2)O2 at Tp

clear; clc;
 
R=8.314*1e-3;   % kJ/molK
Tr=300;         % K
p=1;            % not rellevant
 
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

fprintf('Dissociation temperature: %.2f K \n',Tp);

