function [Tp,n,species,param,flag] = HGSrocket(species,n0,typevec,V0,Pvec,A,options1,options2)
%**************************************************************************
%
% [Tp,n,species,param,flag] = HGSrocket(species,n0,type,V0,P0,P1,A,options)
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
% type --> Entry type that defines the state of the input. 
%          It can be 'T' or 'H'
% V0 --> Entry that should be for type:'T'   V0=T [K] input temperature
%                                      'H'   V0=H [kJ] input enthalpy
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
%                 .maxrange Max range to fit in a parabola
%                 .type Select between: 'Frozen' for frozen flow
%                                       'Shifting' for shifting flow
%               struct('xmin',300,'xmax',5000,'maxiter',50,'epsx',0.1,'epsy',0.5,'fchange',5,'maxrange',1500,
%                   'type','Shifting','info',0,'dTp',100)
%
% Outputs:
%--------------------------------------------------------------------------
% Tp --> [K] Exit temperature
% n --> [mols] Species resultant mols
% species --> String or numbers of species
% param --> Exit parameter: param(1) = Mach
%                           param(2) = Thrust [N]
%                           param(3) = Isp
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

% Rebuild mixtures
if max(id) >= length(HGSdata.Mm)
   [species,n0] = HGSrebuild(species,n0);
   [id] = HGSid(species);
end

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