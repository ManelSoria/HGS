%***********************************************************************************************************
%* HGS 2.0 (adapted from the original HGS 1.3) 
%* Original by Arnau Miro, Pau Manent, Manel Soria 
%* Adapted by Caleb Fuster
%
%* ESEIAAT UPC           
%***********************************************************************************************************
%
% Example 05c: Adiabatic C3H8 reaction using hgsTp
%
% Inlet: H2, O2 
% Outlet: H2O + (1/2)O2 at Tp

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
    [Tp(i),np,~,~]=HGStp(species,nr,'T',298,1,[]);
end

plot(molprop,Tp);
grid;
