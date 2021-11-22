function [Tp,species,n,flag] = HGStp(species,n0,type,V0,P,options)
%**************************************************************************
%
% [Tp,species,n,flag] = HGStp(species,n0,type,V0,P,options)
%
%**************************************************************************
% 
% HGStp calculates the reaction temperature considering 
% dissociation, and the products composition in equilibrium
%
%**************************************************************************
% Inputs:
%--------------------------------------------------------------------------
% species --> String or code of species 
% n0 --> [mols] Number of mols of each species
%        usually, the number of mol of the products will be zero
% type --> Entry type that defines the state of the input. 
%          It can be 'T' or 'H'
% V0 --> Value of type:'T'   V0=T [K] input temperature
%                      'H'   V0=H [kJ] input enthalpy
% P --> [bar] Pressure
% options --> (OPTIONAL) Structure with the options for the secant method. 
%                 .xmin [K] Temperature minimum for the solver;
%                 .xmax [K] Temperature maximum for the solver;
%                 .maxiter Max iterations for the solver;
%                 .epsx Diferential T where the solver reachs the solution;
%                 .epsy Diferential S where the solver reachs the solution;
%                 .fchange T difference where secant method is
%                          changed by bisection method;
%                 .info Detailed info == 1; No info == 0.
%                 .maxrange Max range to fit in a parabola
%                 .dTp Improve the velocity with the approximation of
%                 parabola. +- dTp
%           For instance, by default:
%           struct('xmin',300,'xmax',4000,'maxiter',50,'epsx',0.1,'epsy',0.5,
%                   'fchange',5,'maxrange',1500,'info',0,'dTp',100)
%
% Outputs: 
%--------------------------------------------------------------------------
% Tp --> [K] Exit temperature
% n --> [mols] Species resultant mols
% species --> String or code of species
% flag --> Solver error detection: 
%                 1  Solver has reached the solution
%                -1  Solver failed. Maximum iterations
%                -2  Solver failed. Initial sign change not found
%
%**************************************************************************
% Examples:
% HGStp({'H2','O2','H2O','H','O','OH'}, [2 1 0 0 0 0] , 'T', 400, 10)
% HGStp({'H2','O2','H2O','H','O','OH'}, [2 1 0 0 0 0] , 'H', 3.1, 10)
%**************************************************************************
% *HGS 2.1
% *By Caleb Fuster, Manel Soria and Arnau Miró
% *ESEIAAT UPC    

% type = 'H' or 'T'
if type~='H' && type~='T'
    error('Wrong type = %s',type);
end

if ~exist('options','var')
   options = []; 
end

global HGSdata; HGSload
[id] = HGSid(species);

% Rebuild mixtures
if max(id) >= length(HGSdata.Mm)
   [species,n0] = HGSrebuild(species,n0);
   [id] = HGSid(species);
end

% Compute initial enthalpy
if type=='T'
    H = HGSprop(id,n0,V0,P,'H');
else
    H = V0;
end

[Tp,~,n,flag]=HGSeqcond(id,n0,'H',H,P,'Shifting',options);

end
