%***********************************************************************************************************
%* HGS 2.0 (adapted from the original HGS 1.3) 
%* Original by Arnau Miro, Pau Manent, Manel Soria 
%* Adapted by Caleb Fuster
%
%* ESEIAAT UPC              
%***********************************************************************************************************
%
% Example 05b: Adiabatic H2 / O2 reaction using hgsTp
%
% Inlet: H2, O2 
% Outlet: H2O + (1/2)O2 at Tp

clear

species={'H2','O2','H2O','H','O','OH'};
Tr=350; % K 
P=10; % bar
nr=[2;1;0;0;0;0]; % mol

% Solver options
options = [];

[Tp,np,~,~] = HGStp(species,nr,'T',298,1,options)
np/sum(np)