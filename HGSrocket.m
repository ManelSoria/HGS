function [Tp,n,species,param,flag] = HGSrocket(species,n0,type,V0,P0,P1,A,options)
%**************************************************************************
%
% [Tp,n,species,param,flag] = HGSrocket(species,n0,type,V0,P0,P1,A,options)
%
%**************************************************************************
% 
% HGSrocket calculates the temperatures, mixture and exit parameter (v,M
% T and Isp)
%
%**************************************************************************
% Inputs:
%--------------------------------------------------------------------------
% species --> String or numbers of species
% n0 --> [mols] Number of mols of each species
% type --> 
% V0 --> [K] Initial temperature
% P0 --> [bar] Inlet pressure
% P1 --> [bar] Nozzle exit pressure
% A --> [m^2] Nozzle exit area
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
% param --> Exit parameter: param(1) = v [m/s]
%                           param(2) = Mach
%                           param(3) = Thrust [N]
%                           param(4) = Isp
% flag --> Solver error detection: 
%                 1  Solver has reached the solution
%                -1  Solver failed. Maximum iterations
%                -2  Solver failed. Initial sign change not found
%
%**************************************************************************
% *HGS 2.0
% *By Caleb Fuster, Manel Soria and Arnau Miró
% *ESEIAAT UPC    

% Combustion chamber
[Tp(1),n(:,1),species,flag(1)] = HGStp(species,n0,type,V0,P,options);

% Nozzle
if flag(1) == 1
    [Tp(2),n(:,2),~,param(1),param(2),flag(2)] = HGSisentropic(species,n0,T0,P0,P1,options);
else
    return
end

% Parameters
[MM] = HGSprop(species,n,[],[],'Mm');
m = MM*sum(n);
param(3) = m*param(1)-A*(P1-Pa);
g0 = 9.807; % [m/s^2]
param(4) = param(1)/g0;

end