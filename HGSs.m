function [S] = HGSs(a,T,P,state)
%**************************************************************************
%
% [S] = HGSs(a,T,P,state)
%
%**************************************************************************
%
% HGSh calculates the enthropy of a species using his Burcat coeficients
% and temperature
%
%**************************************************************************
% Inputs:
%--------------------------------------------------------------------------
% a -->  Burcat coefficients
% T --> [K] Temperature
% P --> [bar] Pressure
% state -->  State of the species
%
% Outputs:
%--------------------------------------------------------------------------
% S --> [kJ/(mol*K)]Entropy
%
%**************************************************************************
% *HGS 2.0
% *By Caleb Fuster, Manel Soria and Arnau Miró
% *ESEIAAT UPC    

global R; HGSr

Pref = 1; % [bar]

S = R * (a(7) + a(1)*log(T) + sum(a(2:5).*(T.^(1:4))./(1:4))); % [kJ/mol*K]

if strcmp('G',state) &&  P ~= 0
    S = S-R*log(P/Pref); 
end

end