function [g] = HGSg(s,h,T)
%**************************************************************************
%
% [g] = HGSg(s,h,T)
%
%**************************************************************************
%
% HGSg calculates the species free Gibbs energy using enthalpy, enthropy 
% and temperature
%
%**************************************************************************
% Inputs:
%--------------------------------------------------------------------------
% s --> [kJ/(mol*K)] Entropy
% h --> [kJ/mol] Enthalpy
% T --> [K] Temperature
%
% Outputs:
%--------------------------------------------------------------------------
% g --> [kJ/mol] Free Gibbs energy
%
%**************************************************************************
% *HGS 2.1
% *By Caleb Fuster, Manel Soria and Arnau Miró
% *ESEIAAT UPC    

g = h - T*s; % [kJ/mol]

end