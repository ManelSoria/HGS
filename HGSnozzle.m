function [species,throat,exit,flag] = HGSnozzle(species,n0,T0,P0,Pe,Pa,Fro_Shift,options1,options2)
%**************************************************************************
%
% [species,throat,exit,flag] = HGSnozzle(species,n0,T0,P0,Pe,Pa,Fro_Shift,
%                                           options1,options2)
%
%**************************************************************************
% 
% HGSnozzle calculates the throat and exit conditions from the inlet nozzle
%  conditions
%  
%**************************************************************************
% Inputs:
%--------------------------------------------------------------------------
% species --> String or code of inlet species
% n0 --> [mols] Number of mols/s of each inlet species
% T0 --> [K] Inlet temperature
% P0 --> [bar] Inlet pressure
% Pe --> [bar] Exit pressure
% Pa --> [bar] Atmospheric pressure
% Fro_Shift --> Select between: 'Frozen' for frozen flow
%                               'Shifting' for shifting flow
%                               'Combined' for shifting flow until throat
%                                       and frozen from throat to exit 
% options1 --> Structure with the options for the secant method. (optional)
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
% options2 --> Structure with the options as options1 but for Pressure (optional)
%           struct('xmin',0.01,'xmax',<P0,'maxiter',50,'epsx',0.01,'epsy',0.01,
%                   'fchange',1,'info',0)
%                   Note xmax must be equal or less than inlet Pressure (P0)
%
% Outputs:
%--------------------------------------------------------------------------
% species --> String or code of species
% throat --> cell contaning
%           {1} nt --> [mols] Species mols in throat
%           {2} Tt --> [K] Throat temperature
%           {3} Pt --> [bar] Throat pressure
%           {4} At --> [m^2] Throat area
% exit --> cell contaning
%           {1} ne --> [mols] Species mols in exit
%           {2} Te --> [K] Exit temperature
%           {3} Me --> [] Exit Mch
%           {4} Ae --> [m^2] Exit area
%           {5} F --> [N] Thrust
%           {6} Isp --> [s^]Specific impulse, g0 = 9.807 m/s^2
% flag --> Solver error detection: 
%                 1  Solver has reached the solution
%                -1  Solver failed. Maximum iterations
%                -2  Solver failed. Initial sign change not found
%                Only for Fro_shift = 'Combined'
%                -3  Solver failed. Maximum iterations in P loop
%                -4  Solver failed. Initial sign change not found in P loop
%
%**************************************************************************
% *HGS 2.1
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

% Check input 
if T0<800 || P0<10 
    error('HGSnozzle T0=%e P0=%e too low inlet pressure or temperature',T0,P0);
end

% Total mass
[mm] = HGSprop(id,n0,T0,P0,'Mm'); % g/mol
m=sum(n0)*mm*1e-3; % kg/s


% Throat calculation
if strcmp(Fro_Shift,'Shifting') || strcmp(Fro_Shift,'Combined')
    [Tt,~,nt,Pt,flag] = HGSisentropic(id,n0,T0,P0,'Shifting','M',1,options1,options2);
    if flag ~=1, error('HGSnozzle failed to converge/1 flag=%d',flag), end        
elseif  strcmp(Fro_Shift,'Frozen')
    [Tt,~,nt,Pt,flag] = HGSisentropic(id,n0,T0,P0,'Frozen','M',1,options1,options2);
    if flag ~=1, error('HGSnozzle failed to converge/2 flag=%d',flag), end    
else
    error('Your variable Fro_Shift is no one accepted by this function. Only Frozen and Shifting are accepted')
end

% Exit calculation
if strcmp(Fro_Shift,'Shifting') || strcmp(Fro_Shift,'Frozen')
    % Shifting and Frozen doesnt require extra action
    [Te,~,ne,Me,flag] = HGSisentropic(id,n0,T0,P0,Fro_Shift,'P',Pe,options1);   
    if flag ~=1, error('HGSnozzle failed to converge/3 flag=%d',flag), end    
else
    % Going back until the nozzle start but with a frozen composition 
    [v2,H2] = HGSprop(id,nt,Tt,Pt,'a','H');
    H1 = H2 +(m/2000)*v2^2; % Enthalpy of the start v=0
    [TinitF,~,~,flag]=HGSeqcond(id,nt,'H',H1,P0,'Frozen',options1);
    if flag ~=1, error('HGSnozzle failed to converge/4 flag=%d',flag), end    

    % Frozen expansion
    [Te,~,ne,Me,flag] = HGSisentropic(id,nt,TinitF,P0,'Frozen','P',Pe,options1);
    if flag ~=1, error('HGSnozzle failed to converge/5 flag=%d',flag), end    
end



% Other properties of the nozzle throat
[Rg,a,Mm] = HGSprop(id,nt,Tt,Pt,'Rg','a','Mm');
% Ideal Gas law
rho = Pt/(Rg*Tt);
% Continuity
At = sum(nt)/(Mm*1e-3*a*rho);


% Other properties of the nozzle exit
[Rg,a,Mm] = HGSprop(id,nt,Te,Pe,'Rg','a','Mm');
v2=Me*a;
% Ideal Gas law
rho = Pt/(Rg*Tt);
% Continuity
Ae = sum(ne)/(Mm*1e-3*v2*rho); % Area of the nozzle (m^2)

% Properties
F = m*v2-Ae*(Pe-Pa);
g0 = 9.807;
Isp = v2/g0;

throat = {nt, Tt, Pt, At};
exit = {ne, Te, Me, Ae, F, Isp};




end