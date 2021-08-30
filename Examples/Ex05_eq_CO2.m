%***********************************************************************************************************
% *HGS 2.0
% *By Caleb Fuster, Manel Soria and Arnau Miró
% *ESEIAAT UPC              
%***********************************************************************************************************
%
%  CO2 dissociation equilibrium solved with K
%
% 1 CO2    <-> 1 CO + (1/2)O2 
%   1-z          z       z/2

clear; 

Ru=8.3144621*1e-3;   % kJ/molK
T=2500;         % K
p=1;           % bar


pref=1; % bar (arbitrary pressure)

% Gibbs free energy (g) is computed for each species separately. The rest
% of properties are ignored.
gCO2 = HGSprop({'CO2'},1,T,pref,'G'); % kJ
gCO = HGSprop({'CO'} ,1,T,pref,'G');
gO2 = HGSprop({'O2'} ,1,T,pref,'G');

% Definition of K
deltag=gCO+0.5*gO2-gCO2;
K=exp(-deltag/(Ru*T));

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

% The equation to solve is
eq=@(z) lhs(z)-K;

% With the following parameters of the solver
options=optimset(...
        'Display','none',...
        'MaxIter',4000,...
        'TolFun', 1.0e-10,...
        'TolX',1.0e-4);
 
zi=0.5; % a initial value (arbitrary) of z is given
[z,fval,exitflag]=fsolve(eq,zi,options);

fprintf('K=%e z=%e\n',K,z);

fprintf('x_CO =%6.4f \n',xCO(z));
fprintf('x_CO2=%6.4f \n',xCO2(z));
fprintf('x_O2 =%6.4f \n',xO2(z));