function [Tp,species,n,flag] = HGSeqcond(species,n0,type,V0,P,Fro_Shift,options)
%**************************************************************************
%
% [Tp,species,n,flag] = HGSeqcond(species,n0,type,V0,P,options)
%
%**************************************************************************
% 
% HGSeqcond computes temperature and reaction products that satisfy one of
% the following conditions
% if type=='H', the temperature Tp satisfies the condition H(Tp) = V0
% if type=='S', the temperature Tp satisfies the condition S(Tp) = V0
%
%**************************************************************************
% Inputs:
%--------------------------------------------------------------------------
% species --> Species name, ids are also accpeted as entry
% n0 --> [mols] Number of mols of each species
% type --> Entry type. It coould be 'H' or 'S'
% V0 --> Entry that should be for type:'H'   V0=H [kJ]
%                                      'S'   V0=S [kJ/K]
% P --> [bar] Mixture pressure
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
%                 .maxrange Max range to fit in a parabola
%                 .dTp Improve the velocity with the approximation of
%                 parabola. +- dTp
%           struct('xmin',300,'xmax',6000,'maxiter',50,'epsx',0.1,'epsy',0.5,'fchange',5,'maxrange',1500,
%                   'type','Shifting','info',0,'dTp',100)
%
% Outputs:
%--------------------------------------------------------------------------
% Tp --> [K] Exit temperature
% n --> [mols] Species resultant mols
% species --> String or numbers of species
% flag --> Solver error detection: 
%                 1  Solver has reached the solution
%                -1  Solver failed. Maximum iterations
%                -2  Solver failed. Initial sign change not found
%
% *************************************************************************
% *HGS 2.1
% *By Caleb Fuster, Manel Soria and Arnau Miró
% *ESEIAAT UPC    


if type~='H' && type~='S'
    error('Wrong type = %s',type);
end

global HGSdata; HGSload
[id] = HGSid(species);

% Rebuild mixtures
if max(id) >= length(HGSdata.Mm)
   [species,n0] = HGSrebuild(species,n0);
   [id] = HGSid(species);
end

%% Options

def.xmin = 300;
def.xmax = 4000;
def.maxiter = 200;
def.epsx = 0.1;
def.epsy = 1;
def.fchange = 5;
def.info = 0;
def.dTp = 100;

if exist('options','var') && ~isempty(options)
    fields = fieldnames(options);
    for ii = 1:length(fields)
        def.(fields{ii}) = options.(fields{ii});
    end
end

options = def;

if ~exist('Fro_Shift','var') || strcmpi(Fro_Shift,'')
    Fro_Shift = 'Shifting';
end

%% Secant method resolution

[Tp,n,flag] = HGSsecant(@hastobezero,n0,options);


    function [H,n] = hastobezero(T,n)
        if strcmpi(Fro_Shift,'Shifting')
            [~,n,~]=HGSeq(id,n,T,P,[]);
        end
        [H] = HGSprop(id,n,T,P,type)-V0;
    end

    

end


