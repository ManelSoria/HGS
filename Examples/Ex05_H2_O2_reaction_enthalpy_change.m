%***********************************************************************************************************
%* HGS 2.0 (adapted from the original HGS 1.3) 
%* Original by Arnau Miro, Pau Manent, Manel Soria 
%* Adapted by Caleb Fuster
%
%* ESEIAAT UPC    
%***********************************************************************************************************
%
% Example 05: Adiabatic H2 / O2 reaction
%
% Inlet: H2, O2 
% Outlet: H2O + (1/2)O2 at Tp

clear
close 

species={'H2','O2','H2O','H','O','OH'};
Tr=350; % K 
P=66; % bar
nr=[2;1;0;0;0;0]; % mol
 
 Hin = HGSprop(species,nr,Tr,P,'H')
 
T=linspace(300,5000,30);
for i=1:length(T) 
    [~,comp,~]= HGSeq(species,nr,T(i),P);
    Hout(i) = HGSprop(species,comp,T(i),P,'H');    
end
 
plot(T,Hout,'-or',T,Hin*ones(length(T),1),'-ob')