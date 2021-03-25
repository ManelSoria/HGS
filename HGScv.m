function [Cv] = HGScv(Cp)
%**************************************************************************
%
% [Cv] = HGScv(Cp)
%
%**************************************************************************
%
% HGScv calculates the species Constant volume coeficient using constant
% pressure coeficient
%
%**************************************************************************
%
% Inputs:
%--------------------------------------------------------------------------
% Cp --> [kJ/(mol*K)] Constant pressure coeficient
% Mm --> [g/mol] Molar Mass
%
% Outputs:
%--------------------------------------------------------------------------
% Cv --> [kJ/(mol*K)] Constant volume coeficient
%
%**************************************************************************
% *HGS 2.0
% *By Caleb Fuster, Manel Soria and Arnau Miró
% *ESEIAAT UPC    

global R

Cv = Cp - R;  

end