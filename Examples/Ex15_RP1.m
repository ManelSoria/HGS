%***********************************************************************************************************
% *HGS 2.0
% *By Caleb Fuster, Manel Soria and Arnau Miró
% *ESEIAAT UPC         
%***********************************************************************************************************
%
% Operating with mixtures: RP1 surrogate model and air

clear; 

format compact
% lets see what are RP1 and Air7921
[RP1species,RP1n]=HGSrebuild('RP1',1)
[Air7921,Air7921n]=HGSrebuild('Air7921',1)


% now RP1 can be used as a single species
species={'RP1','O2','CO2','CO','H2O','OH','O','H'};
n      =[1    , 1  , 0   , 0  , 0   , 0  , 0 , 0];
[Mm,Cp,Cv,H,S,G,Rg,gamma,a]=HGSprop(species,n,1500,100); 

% RP1 and air combustion[
[T,n,sp,flag]=HGStp({'RP1','Air7921','CO2','CO','H2O','OH','O','H'},[1, 1, 0, 0, 0, 0, 0, 0 ],'T',300,40) % not O2 as it is already in Air7921 

