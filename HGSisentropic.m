function [Tp,species,n,V2,flag] = HGSisentropic(species,n0,T0,P0,Fro_Shift,typeexit,V1,options1,options2)
%**************************************************************************
%
% [Tp,species,n,V2,flag] = HGSisentropic(species,n0,T0,P0,Fro_Shift,
%                                        typeexit,V1,options1,options2)
%
%**************************************************************************
% 
% HGSisentropic calculates the outlet variables for an isentropic expansion 
%
%**************************************************************************
% Inputs:
%--------------------------------------------------------------------------
% species --> String or code of species
% n0 --> [mols] Number of mols of each species
% T0 --> [K] Initial temperature
% P0 --> [bar] Inlet pressure
% Fro_Shift --> Select between: 'Frozen' for frozen flow
%                               'Shifting' for shifting flow
% typeexit --> Entry type that defines the state of the input. 
%               It can be 'P' or 'M'
% V1 --> Value for type:'P'   V1=P [bar] output pressure
%                       'M'   V1=M [] output Mach. Has to be >=1
% options1 --> (optional) Structure with the options to be passed to HGSeqcond 
%                  that solves for the temperature with HGSsecant. 
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
%           struct('xmin',300,'xmax',4000,'maxiter',50,'epsx',0.1,'epsy',0.5,
%                   'fchange',5,'maxrange',1500,'info',0,'dTp',100)
% options2 --> (optional) Structure with the options for HGSsecant to be 
%              called here to solve for the Pressure
%              struct('xmin',0.01,'xmax',<P0,'maxiter',50,'epsx',0.01,'epsy',0.01,
%                   'fchange',1,'info',0)
%
%**************************************************************************
% Outputs:
%--------------------------------------------------------------------------
% Tp --> [K] Exit temperature
% species --> String or code of species
% n --> [mols] Species c mols
% V2 --> If typeexit 'P' ->[Mach] Mach of the mixture at a P
%        If typeexit 'M' ->[bar] Pressure at the Mach
% flag --> Solver error detection: 
%                 1  Solver has reached the solution
%                -1  Solver failed. Maximum iterations in T loop
%                -2  Solver failed. Initial sign change not found in T loop
%                Only for typeexit=='M'
%                -3  Solver failed. Maximum iterations in P loop
%                -4  Solver failed. Initial sign change not found in P loop
%
%**************************************************************************
% *HGS 2.1
% *By Caleb Fuster, Manel Soria and Arnau Miró
% *ESEIAAT UPC    
  

global HGSdata; HGSload
[id] = HGSid(species);

if ~exist('options1','var')
   options1 = struct('xmin',300,'xmax',5000,'maxiter',50,'epsx',1,'epsy',10^(ceil(log10(max(n0)))-1)*1e-4,'fchange',25,'maxrange',1500,'info',0,'dTp',100); 
end

% Rebuild mixtures
if max(id) >= length(HGSdata.Mm)
   [species,n0] = HGSrebuild(species,n0);
   [id] = HGSid(species);
end

% Compute initial entropy and enthalpy, mm

[S,Mm1,H1] = HGSprop(id,n0,T0,P0,'S','Mm','H');% Inlet properties
m=sum(n0)*Mm1*1e-3;
h1=H1/m;

if strcmpi(typeexit,'P')
    P1 = V1;
    [Tp,~,n,flag]=HGSeqcond(id,n0,'S',S,V1,Fro_Shift,options1); % Searching T so that S=S0 and P=P1
    [a2,H2] = HGSprop(id,n,Tp,P1,'a','H'); % Outlet properties

    h2=H2/m;

    v2=sqrt(2*1000*(h1-h2)); % Enthalpy must be en J/kg !

    V2=v2/a2;
elseif strcmpi(typeexit,'M')
    % To avoid problems with imaginary numbers
    if strcmpi('Frozen',Fro_Shift)
        P0=P0*0.9; 
    end
    % Change options to solve for P, imposing M
    if ~exist('options2','var') || isempty(options2)
        options2 = struct('xmin',0.1,'xmax',P0,'maxiter',50,'epsx',0.01,'epsy',0.001,'fchange',1,'info',0);
    else
        options2 = UpdateOpt(options2);
    end
    % Pressure search, imposing M
    [V2,n,flag] = HGSsecant(@hastobeM,n0,options2);
    % Flag error return
    if flag~=1
        Tp=[];n=[];V2=[];flag = flag-2;
        return
    end    
    % T calculation
    [Tp,~,n,flag]=HGSeqcond(id,n,'S',S,V2,Fro_Shift,options1);
else
    error('Typeexit is not accepted by this function. Only accepted P and M')
end

    function [zeroM,n] = hastobeM(Pstar,n)
        function [zeroS,n] = hastobeS(Tstar,n)
            if strcmpi(Fro_Shift,'Shifting')
               [species,n,~] = HGSeq(species,n,Tstar,Pstar);
            end
            S1 = HGSprop(species,n,Tstar,Pstar,'S'); 
            zeroS = (S1 - S);   
        end          
        [Tstar,n,flags] = HGSsecant(@hastobeS,n,options1);
        if flags~=1
            error('uhhh error in hastobeM HGSsecant failed flag=%d \n',flags);
        end
        [a,H2] = HGSprop(species,n,Tstar,Pstar,'a','H');   
        v2=sqrt(2*1000*(h1-H2/m));
        M1 = v2/a;
        zeroM = V1 - M1;
    end
    

    function opt = UpdateOpt(options)
        %Overwrite default options with the new ones
        opt.xmin = 0.1;
        opt.xmax = P0;
        opt.maxiter = 50;
        opt.epsx = 0.01;
        opt.epsy = 0.001;
        opt.fchange = 50;
        opt.info = 1;

        if ~isempty(options)
            fields = fieldnames(options);
            for ii = 1:length(fields)
                opt.(fields{ii}) = options.(fields{ii});
            end
        end

    end

end