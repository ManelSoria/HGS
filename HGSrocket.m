function [Tp,n,species,param,flag] = HGSrocket(species,n0,typevec,V0,Pvec,A,options1,options2)
%**************************************************************************
%
% [Tp,n,species,param,flag] = HGSrocket(species,n0,typevec,V0,Pvec,A,
%                                       options1,options2)
%
%**************************************************************************
% 
% HGSrocket calculates the temperatures, mixture and exist parameter (v,M
% T and Isp)
%
%**************************************************************************
% Inputs:
%--------------------------------------------------------------------------
% species --> String or numbers of species
% n0 --> [mols] Number of mols of each species
% typevec --> {1} Entry type that defines the state of the combustion input. 
%                 It can be 'T' or 'H'
%             {2} Select between: 'Frozen' for frozen flow
%                                 'Shifting' for shifting flow
%                                 'Combined' for shifting flow until throat
%                                           and frozen until end 
% V0 --> Entry that should be for typevec{1}:'T'   V0=T [K] input temperature
%                                            'H'   V0=H [kJ] input enthalpy
% Pvec --> (1)[bar] Inlet pressure
%          (2)[bar] Nozzle exit pressure
%          (3)[bar] Atmospheric pressure
% A --> [m^2] Nozzle exit area
% options --> Structure with the options for the secant method. 
%                 .xmin [K] Temperature minimum for the solver;
%                 .xmax [K] Temperature maximum for the solver;
%                 .maxiter Max iterations for the solver;
%                 .epsx Diferential T where the solver reachs the solution;
%                 .epsy Diferential S where the solver reachs the solution;
%                 .fchange T difference where secant method is
%                          changed by bisection method;
%                 .maxrange Max range to fit in a parabola
%                 .type Select between: 'Frozen' for frozen flow
%                                       'Shifting' for shifting flow
%               struct('xmin',300,'xmax',5000,'maxiter',50,'epsx',0.1,'epsy',0.5,'fchange',5,'maxrange',1500,
%                   'type','Shifting','info',0,'dTp',100)
% options2 --> Structure with the options as options1 but for Pressure
%           struct('xmin',0.01,'xmax',<P0,'maxiter',50,'epsx',0.01,'epsy',0.01,
%                   'fchange',1,'info',0)
%
% Outputs:
%--------------------------------------------------------------------------
% Tp --> (1)[K] Exit temperature of the CC
%        (2)[K] Exit temperature of the nozzle
% n --> (1)[mols] Species resultant mols of the CC
%       (2)[mols] Species resultant mols of the nozzle    
% species --> String or numbers of species
% param --> Exit parameter: param(1) = Mach
%                           param(2) = Thrust [N]
%                           param(3) = Isp
% flag --> Solver error detection: 
%                 1  Solver has reached the solution
%                -1  Solver failed. Maximum iterations
%                -2  Solver failed. Initial sign change not found
%                Only for typevec{2}=='Combinated'
%                -3  Solver failed. Maximum iterations in P loop
%                -4  Solver failed. Initial sign change not found in P loop
%
%**************************************************************************
% *HGS 2.0
% *By Caleb Fuster, Manel Soria and Arnau Miró
% *ESEIAAT UPC    

global HGSdata; HGSload
[id] = HGSid(species);

% Rebuild mixtures
if max(id) >= length(HGSdata.Mm)
   [species,n0] = HGSrebuild(species,n0);
   [id] = HGSid(species);
end

% Options
if ~exist('options1','var')
    options1 = [];
end
if ~exist('options2','var')
    options2 = [];
end


% Combustion chamber
[Tp(1),n(1,:),species,flag(1)] = HGStp(id,n0,typevec{1},V0,Pvec(1),options1);

if flag(1)~=1
   return 
end

% Nozzle
[Tp(2),n(2,:),~,param(1),param(2),param(3),flag(2)] = HGSnozzle(id,n(1,:),Tp(1),Pvec(1),Pvec(2),Pvec(3),A,typevec{2},options1,options2);

end