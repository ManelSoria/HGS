%***********************************************************************************************************
% *HGS 2.0
% *By Caleb Fuster, Manel Soria and Arnau Miró
% *ESEIAAT UPC          
%***********************************************************************************************************
%
% Adiabatic C3H8 reaction using hgsTp
%
% Inlet: H2, O2 
% Outlet: H2O + (1/2)O2 at Tp
% TBD: rewrite the example so that the variable is:
% (1) the OF ratio 
% (2) the O2 excess

clear
species={'C3H8',...
    'CO2',...
    'CO',...
    'O2',...
    'O',...
    'H2',...
    'H',...
    'OH',...
    'H2O'};

Tr=280; % K 
P=5; % bar

molprop=linspace(1,4,15);
for i=1:length(molprop)
    fprintf('Solving for i=%d / %d ... \n',i,length(molprop)); 
    nr=[molprop(i);0;0;5;0;0;0;0;0]; % mol
    [Tp(i),~,np,~]=HGStp(species,nr,'T',298,1);
end

plot(molprop,Tp,'LineWidth',2);
grid;
set(gca,'FontSize',18)
xlabel('n (mol C3H8)');
ylabel('T (K)');
title('Flame temperature (with dissociation) of n mol C3H8 + 5 mol O2');

