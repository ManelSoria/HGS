%***********************************************************************************************************
% *HGS 2.1
% *By Caleb Fuster, Manel Soria and Arnau Miró
% *ESEIAAT UPC 
%***********************************************************************************************************
%
% H2O equilibrium dissociation for different values of
%             temperature (T)
%
% H20 <-> H2 + O2 + H + O + OH

clear; clc;

p=1;                        % bar
T=linspace(300,5000,50);    % K

% loop to compute composition using hgseq
lenT = length(T);
xH2 = zeros(lenT,1);
xO2 = zeros(lenT,1);
xH2O = zeros(lenT,1);
xH = zeros(lenT,1);
xO = zeros(lenT,1);
xOH = zeros(lenT,1);
for i=1:lenT
    fprintf('Solving equilibrium composition for T=%f K\n',T(i));
    [~,comp,~]=HGSeq({'H2','O2','H2O','H','O','OH'},[2;1;0;0;0;0],T(i),1);
    xH2(i)=comp(1)/sum(comp);
    xO2(i)=comp(2)/sum(comp);
    xH2O(i)=comp(3)/sum(comp);
    xH(i)=comp(4)/sum(comp);
    xO(i)=comp(5)/sum(comp);
    xOH(i)=comp(6)/sum(comp);
end

% plot the results and apply some basic formatting
plot(T,xH2O,'r',T,xH,'b',T,xO,'g',T,xOH,'k',T,xH2,'c',T,xO2,'m','LineWidth',2);
legend('H20','H','O','OH','H2','O2','Location','NorthWest');
xlabel('Temperature (K)'); ylabel('Molar fraction');
set(gca,'FontSize',18)
title('H_2O dissociation as a function of temperature');
grid;
