function [G] = HGSg(S,H,T)
%**************************************************************************
%
% [G] = HGSg(S,H,T)
%
%**************************************************************************
%
% HGSg calculates the species free Gibbs energy using enthalpy, enthropy 
% and temperature
%
%**************************************************************************
% Inputs:
%--------------------------------------------------------------------------
% S --> [kJ/(mol*K)] Enthropy
% H --> [kJ/mol] Enthalpy
% T --> [K] Temperature
%
% Outputs:
%--------------------------------------------------------------------------
% G --> [kJ/mol] Free Gibbs energy
%
%**************************************************************************
% *HGS 2.0
% *By Caleb Fuster, Manel Soria and Arnau Miró
% *ESEIAAT UPC    

G = H - T*S; % [kJ/mol]

end