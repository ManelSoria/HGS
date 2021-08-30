function [Tp,n,species,v2,M2,flag] = HGSisentropic(species,n0,T0,P0,P1,options)
%**************************************************************************
%
% [Tp,n,species,v2,M2,flag] = HGSisentropic(species,n0,T0,P0,P1,options)
%
%**************************************************************************
% 
% HGSisentropic calculates the outlet variables for an isentropic expansion 
%
%**************************************************************************
% Inputs:
%--------------------------------------------------------------------------
% species --> String or numbers of species
% n0 --> [mols] Number of mols of each species
% T0 --> [K] Initial temperature
% P0 --> [bar] Inlet pressure
% P1 --> [bar] Exit pressure
% options --> Structure with the options for the secant method. 
%                 .xmin [K] Temperature minimum for the solver;
%                 .xmax [K] Temperature maximum for the solver;
%                 .maxiter Max iterations for the solver;
%                 .epsx Diferential T where the solver reachs the solution;
%                 .epsy Diferential S where the solver reachs the solution;
%                 .fchange T difference where secant method is
%                          changed by bisection method;
%                 .type Select between: 'Frozen' for frozen flow
%                                       'Shifting' for shifting flow
%                 .info Detailed info == 1; No info == 0.
%                 .dTp Improve the velocity with the approximation of
%                 parabola. +- dTp
%           struct('xmin',300,'xmax',6000,'maxiter',50,'epsx',0.1,'epsy',0.5,'fchange',5,'type','Shifting','info',0,'dTp',100)
%
%**************************************************************************
% Outputs:
%--------------------------------------------------------------------------
% Tp --> [K] Exit temperature
% n --> [mols] Species resultant mols
% species --> String or numbers of species
% v2 --> [m/s] Velocity of the mixture
% M2 --> [Mach] Mach of the mixture
% flag --> Solver error detection: 
%                 1  Solver has reached the solution
%                -1  Solver failed. Maximum iterations
%                -2  Solver failed. Initial sign change not found
%
%**************************************************************************
% *HGS 2.0
% *By Caleb Fuster, Manel Soria and Arnau Miró
% *ESEIAAT UPC    
  

global HGSdata; HGSload
[id] = HGSid(species);

if ~exist('options','var')
   options = []; 
end

% Rebuild mixtures
if max(id) >= length(HGSdata.Mm)
   [species,n0] = HGSrebuild(species,n0);
   [id] = HGSid(species);
end

% compute initial entropy and enthalpy

[S,Mm1,H1] = HGSprop(id,n0,T0,P0,'S','Mm','H');% Inlet properties

[Tp,n,~,flag]=HGSeqcond(id,n0,'S',S,P1,options);

[Mm2,a2,H2] = HGSprop(id,n,Tp,P1,'Mm','a','H'); % Outlet properties

m1=sum(n0)*Mm1*1e-3;
h1=H1/m1;
m2=sum(n)*Mm2*1e-3;
h2=H2/m2;

v2=sqrt(2*1000*(h1-h2)); % Enthalpy must be en J/kg !

M2=v2/a2;

end