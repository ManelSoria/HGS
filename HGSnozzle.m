function [Tp,np,species,M,F,Isp,flag] = HGSnozzle(species,n0,T0,P0,P1,Pa,A,Fro_Shift,options1,options2)
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
% Pa --> [bar] Atmospheric pressure
% A --> [m^2] Exit area
% Fro_Shift --> Select between: 'Frozen' for frozen flow
%                               'Shifting' for shifting flow
%                               'Combined' for shifting flow until throat
%                                       and frozen until end 
% options1 --> Structure with the options for the secant method. 
%                 .xmin [K] Temperature minimum for the solver;
%                 .xmax [K] Temperature maximum for the solver;
%                 .maxiter Max iterations for the solver;
%                 .epsx Diferential T where the solver reachs the solution;
%                 .epsy Diferential S where the solver reachs the solution;
%                 .fchange T difference where secant method is
%                          changed by bisection method;
%                 .maxrange Max range to fit in a parabola
%                 .info Detailed info == 1; No info == 0.
%                 .dTp Improve the velocity with the approximation of
%                 parabola. +- dTp
%           struct('xmin',300,'xmax',6000,'maxiter',50,'epsx',0.1,'epsy',0.5,
%                   'fchange',5,'maxrange',1500,'info',0,'dTp',100)
% options2 --> Structure with the options as options1 but for Pressure
%           struct('xmin',0.01,'xmax',<P0,'maxiter',50,'epsx',0.01,'epsy',0.01,
%                   'fchange',1,'info',0)
%
% Outputs:
%--------------------------------------------------------------------------
% Tp --> [K] Exit temperature
% np --> [mols] Species resultant mols
% species --> String or numbers of species
% M --> [] Exit Mach
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

if ~exist('options1','var')
    options1 = [];
end
if ~exist('options2','var')
    options2 = struct('xmin',P1,'xmax',P0,'maxiter',50,'epsx',0.01,'epsy',0.001,'fchange',1,'info',0);
end

[mm] = HGSprop(species,n0,T0,P0,'Mm');
m=sum(n0)*mm*1e-3;

if strcmp(Fro_Shift,'Shifting') || strcmp(Fro_Shift,'Frozen')
    [Tp,np,species,M,flag] = HGSisentropic(species,n0,T0,P0,Fro_Shift,'P',P1,options1);
elseif strcmp(Fro_Shift,'Combinated')
    [Tt,nt,species,Pt,flag] = HGSisentropic(species,n0,T0,P0,'Shifting','M',1,options1,options2);
    if flag ~=1
       return 
    end
    [v2,H2] = HGSprop(species,nt,Tt,Pt,'a','H');
    H1 = H2 +(m/2000)*v2^2;
    [TinitF,~,~,flag]=HGSeqcond(species,nt,'H',H1,P0,'Frozen',options1);
    if flag ~=1
       return 
    end
    [Tp,np,species,M,flag] = HGSisentropic(species,nt,TinitF,P0,'Frozen','P',P1,options1);
else
    error('Your variable Fro_Shift is no one accepted by this function. Only Frozen and Shifting are accepted')
end
if flag ~=1
   return 
end

[a] = HGSprop(species,np,Tp,P1,'a');
v2=M*a;
F = m*v2-A*(P1-Pa);
g0 = 9.807;
Isp = v2/g0;


end