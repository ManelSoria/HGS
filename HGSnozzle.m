function [Tp,n,species,F,Isp,flag] = HGSnozzle(species,n0,T0,P0,P1,A,options)
%**************************************************************************
%
% [Tp,n,species,F,Isp,flag] = HGSnozzle(species,n0,T0,P0,P1,A,options)
%
%**************************************************************************
% 
% HGSnozzle calculates from the inlet nozzle conditions the outlet 
% conditions plus the thrust and Isp. 
%
%**************************************************************************
% Inputs:
%--------------------------------------------------------------------------
% species --> String or numbers of species
% n0 --> [mols] Number of mols of each species
% T0 --> [K] Initial temperature
% P0 --> [bar] Inlet pressure
% P1 --> [bar] Exit pressure
% A --> [m^2] Exit area
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
%           struct('xmin',300,'xmax',6000,'maxiter',50,'epsx',0.1,'epsy',0.5,'fchange',5,'type','Shifting','info',0,'')
%
% Outputs:
%--------------------------------------------------------------------------
% Tp --> [K] Exit temperature
% n --> [mols] Species resultant mols
% species --> String or numbers of species
% F --> [N] Thrust
% Isp --> [s^-1]Specific impulse, g0 = 9.807 m/s^2
% flag --> Solver error detection: 
%                 1  Solver has reached the solution
%                -1  Solver failed. Maximum iterations
%                -2  Solver failed. Initial sign change not found
%
%**************************************************************************
% *HGS 2.0
% *By Caleb Fuster, Manel Soria and Arnau Miró
% *ESEIAAT UPC    

[Tp,n,species,v2,~,flag] = HGSisentropic(species,n0,T0,P0,P1,options);
[MM] = HGSprop(species,n,[],[],'Mm');
m = MM*sum(n);
F = m*v2-A*(P1-Pa);
g0 = 9.807;
Isp = v2/g0;


end