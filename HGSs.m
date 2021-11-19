function [s] = HGSs(a,T,P,state)
%**************************************************************************
%
% [s] = HGSs(a,T,P,state)
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
% s --> [kJ/(mol*K)] Entropy
%
%**************************************************************************
% *HGS 2.1
% *By Caleb Fuster, Manel Soria and Arnau Mir�
% *ESEIAAT UPC    

global R; HGSr

Pref = 1; % [bar]

s = R * (a(7) + a(1)*log(T) + sum(a(2:5).*(T.^(1:4))./(1:4))); % [kJ/mol*K]

if strcmp('G',state) &&  P ~= 0
    s = s-R*log(P/Pref); 
end

end