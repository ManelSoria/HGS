function [s] = HGSs(a,T,p,state)
%**************************************************************************
%
% [s] = HGSs(a,T,p,state)
%
%**************************************************************************
%
% HGSs calculates the enthropy of a species using his Burcat coeficients
% and temperature
%
%**************************************************************************
% Inputs:
%--------------------------------------------------------------------------
% a -->  Burcat coefficients
% T --> [K] Temperature
% p --> [bar] Pressure
% state -->  State of the species
%
% Outputs:
%--------------------------------------------------------------------------
% s --> [kJ/(mol*K)] Entropy
%
%**************************************************************************
% *HGS 2.1
% *By Caleb Fuster, Manel Soria and Arnau Miró
% *ESEIAAT UPC    

global R; HGSr

p_ref = 1; % [bar] Reference pressure

s = R * (a(7) + a(1)*log(T) + sum(a(2:5).*(T.^(1:4))./(1:4))); % [kJ/mol*K]

if strcmp('G',state) &&  p ~= 0
    s = s-R*log(p/p_ref); 
end

end