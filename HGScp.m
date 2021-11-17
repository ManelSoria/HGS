function [Cp] = HGScp(a,T)
%**************************************************************************
%
% [Cp] = HGScp(a,T)
%
%**************************************************************************
%
% HGScp calculates the species Constant pressure coeficient using Burcat
% coeficients  and temperature
%
%**************************************************************************
% Inputs:
%--------------------------------------------------------------------------
% a --> Burcat coefficients
% T --> [K] Temperature
%
% Outputs:
%--------------------------------------------------------------------------
% Cp --> [kJ/(mol*K)] Constant pressure coeficient
%
%**************************************************************************
% *HGS 2.1
% *By Caleb Fuster, Manel Soria and Arnau Miró
% *ESEIAAT UPC    

global R; HGSr;

Cp = R * (a(1) + sum(a(2:5).*(T.^(1:4)))); % [kJ/mol*K]


end