%***********************************************************************************************************
%* HGS 2.0 (adapted from the original HGS 1.3) 
%* Original by Arnau Miro, Pau Manent, Manel Soria 
%* Adapted by Caleb Fuster
%
%* ESEIAAT UPC                
%***********************************************************************************************************
%
% Example 03: CO2 dissociation equilibrium solved with K
% For p=1 bar, plot CO equilibrium fraction as a function of T

clear; 

Ru=8.3144621*1e-3;   % kJ/molK

p=1;

pref=1; % bar (arbitrary pressure)

% Amount of mols as a function of z
nCO=@(z) z;
nO2=@(z) z/2;
nCO2=@(z) 1-z;

% Total amount as a function of z
nT=@(z) nCO(z)+nO2(z)+nCO2(z);

% Molar fractions as a function of z
% Note: this code is probably very slow !
% It would be faster to use a single function for lhs(z,p)
xCO=@(z) nCO(z)/nT(z);
xO2=@(z) nO2(z)/nT(z);
xCO2=@(z) nCO2(z)/nT(z);


% From the definition of K the left hand side (lhs) of the equation is
% written
lhs=@(z) (xCO(z)*xO2(z)^0.5 / xCO2(z) ) * (p/pref)^(1+1/2-1);



% With the following parameters of the solver
options=optimset(...
        'Display','none',...
        'MaxIter',4000,...
        'TolFun', 1.0e-10,...
        'TolX',1.0e-4);
 
zi=0.5; % a initial value (arbitrary) of z is given

TT=linspace(500,3000,20);
tic
for i=1:length(TT)
    
    T=TT(i);
    gCO2 = HGSprop({'CO2'},1,T,pref,'G'); % kJ
    gCO = HGSprop({'CO'} ,1,T,pref,'G');
    gO2 = HGSprop({'O2'} ,1,T,pref,'G');

    % Definition of K
    deltag=gCO+0.5*gO2-gCO2;
    K=exp(-deltag/(Ru*T));
    
    % The equation to solve is
    eq=@(z) lhs(z)-K;
    
    [z,fval,exitflag]=fsolve(eq,zi,options);
    
    zi=z; % initial guess for next temperature
    XX(i)=xCO(z);
end
toc


semilogy(TT,XX,'o-');
grid