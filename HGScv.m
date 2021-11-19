function [cv] = HGScv(cp)
%**************************************************************************
%
% [cv] = HGScv(cp)
%
%**************************************************************************
%
% HGScv calculates the species Constant volume coeficient using constant
% pressure coeficient
%
%**************************************************************************
% Inputs:
%--------------------------------------------------------------------------
% cp --> [kJ/(mol*K)] Constant pressure coefficient
%
% Outputs:
%--------------------------------------------------------------------------
% cv --> [kJ/(mol*K)] Constant volume coefficient
%
%**************************************************************************
% *HGS 2.1
% *By Caleb Fuster, Manel Soria and Arnau Miró
% *ESEIAAT UPC    

global R

cv = cp - R;  

end